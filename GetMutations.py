from collections import namedtuple
from firebrowse import fbget
import numpy as np
import os
import shutil
import subprocess
import tarfile


def GetMutations(FirehosePath, MutsigQ, Disease, Output):
    """Generates variables containing mutation events with significance
    'MutsigQ' or greater. Uses Firebrowse, a tool from the Broad Genome Data
    Analysis Center to download Mutsig2CV results from the Broad Institute
    servers and to extract significant events from these results.
    Automatically cleans up results on completion.

    Parameters
    ----------
    FirehosePath : string
        Path to firehose_get executable.
    MutsigQ : double
        A scalar in the range [0, 1] specifying the Mutsig2CV significance
        threshold to use when filtering somatic mutation events.
    Disease : string
        Dataset code to generate somatic mutation profiles for. Can be obtained
        using firehose_get -c.
    Output : string
        Path to be used for temporary downloading and unzipping Mutsig2CV
        files. Downloads and extracted files will be removed from disk on
        cleanup.

    Returns
    -------
    Mutation : named_tuple
        A named tuple containing the following fields:
        'Symbols' - a numpy array containing gene symbols of Mutsig2CV
                    significant mutation events.
        'Binary' - numpy array of mutations where rows correspond to 'Symbols'
                and columns corresponding to 'Barcodes'. A '1' in position
                (i,j) indicates that gene Symbol[i] is mutated in sample
                Barcode[j].
        'Barcodes' - numpy array containing TCGA barcodes for samples.
        'Type' - data type - 'SomaticMutation' in this case.
        'Release' - release data (string) of Mutsig2CV data from Broad GDAC.

    Notes
    -----
    Return values are returned as a namedtuple.
    """

    # make output folder if exists
    if not os.path.isdir(Output):
        os.mkdir(Output)

    # download Mutsig table, parse and get gene symbols
    Chars = fbget.smg(cohort=Disease, tool="MutSig2CV",
                      format="tsv", q=MutsigQ)
    Mutsig = [row.split("\t") for row in Chars.split("\n")]
    GeneCol = Mutsig[0].index("gene")
    Symbols = [str(Mutsig[i][GeneCol]) for i in range(1, len(Mutsig)-1)]
    Symbols.sort()

    # get name of latest firehose data run
    FH = subprocess.Popen([FirehosePath + "firehose_get -r"],
                          stdout=subprocess.PIPE, shell=True)
    (Runs, err) = FH.communicate()
    Latest = [Run for Run in Runs.split("\n") if Run.startswith("stddata")][-1]

    # fetch mutation data from firehose and move to output
    FH = subprocess.Popen(["cd " + Output + "; " +
                           FirehosePath + "firehose_get -b -tasks " +
                           "Mutation_Packager_Calls.Level_3 stddata latest " +
                           Disease], stdout=subprocess.PIPE, shell=True)
    (Out, err) = FH.communicate()

    # find *.tar.gz files in download directory
    Files = []
    for root, dirs, files in os.walk(Output + Latest):
        for file in files:
            if file.endswith(".tar.gz"):
                Files.append(root + "/" + file)

    # extract maf, manifest files
    MAFZip = [File for File in Files if
              "Mutation_Packager_Calls.Level_3" in File]
    Tar = tarfile.open(MAFZip[0])
    ManifestFile = [member for member in Tar.getmembers() if "MANIFEST.txt"
                    in member.name]
    ManifestFile[0].name = os.path.basename(ManifestFile[0].name)
    Tar.extract(ManifestFile[0], path=Output)
    MAFFiles = [member for member in Tar.getmembers() if
                member.name.find(".maf.txt") != -1]
    MAFFilenames = []
    for member in MAFFiles:
        member.name = os.path.basename(member.name)
        Tar.extract(member, path=Output)
        MAFFilenames.append(member.name)

    # read in manifest and get barcodes that were sequenced
    TextFile = open(Output + "MANIFEST.txt", 'r')
    Contents = np.array([line[:-1].split(' ') for line in TextFile])
    Barcodes = list(Contents[:, 1])
    Barcodes = [Barcode.split('.')[0] for Barcode in Barcodes]

    # initialize matrix for mutations
    Binary = np.zeros((len(Symbols), len(Barcodes)), dtype=np.float)
    for File in MAFFilenames:
        Index = Barcodes.index(File.split('.')[0])
        TextFile = open(Output + File, 'r')
        Contents = [line[:-1].split('\t') for line in TextFile]

    # delete silent mutations
    Variant = Contents[0].index('Variant_Classification')
    Contents = np.array(Contents[1:])
    Silent = (Contents[:, Variant] == 'Silent').nonzero()
    Contents = np.delete(Contents, Silent, axis=0)

    # MAP non-silent maf symbols to 'Symbols' and fill in Mutations
    for i in range(Contents.shape[0]):
        if Contents[i, 0] in Symbols:
            Binary[Symbols.index(Contents[i, 0]), Index] = 1

    # cleanup files
    os.remove(Output + "MANIFEST.txt")
    shutil.rmtree(Output + Latest)
    for File in MAFFilenames:
        os.remove(Output + File)

    # build named tuples for outputs
    MutsigTuple = namedtuple('Mutation', ['Symbols', 'Binary', 'Barcodes',
                                          'Type', 'Release'])
    Mutations = MutsigTuple(Symbols, Binary, Barcodes, 'SomaticMutation',
                            Latest)

    # return outputs
    return Mutations
