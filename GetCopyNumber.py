from collections import namedtuple
from HUGOFilter import HUGOFilter
import numpy as np
import os
import shutil
import subprocess
import tarfile


def GetCopyNumber(FirehosePath, Disease, Output, GisticQ=0.25,
                  FilterGenes=None):
    """Generates variables containing gene-level and arm-level copy number
    levels for events with significance 'GisticQ' or greater. Uses Firebrowse,
    a tool from the Broad Genome Data Analysis Center to download GISTIC
    results from the Broad Institute servers and to extract significant events
    from these results. Automatically cleans up results on completion.

    Parameters
    ----------
    FirehosePath : string
        Path to firehose_get executable.
    Disease : string
        Dataset code to generate copy number profiles for. Can be obtained
        using firehose_get -c.
    Output : string
        Path to be used for temporary downloading and unzipping GISTIC files.
        Downloads and extracted files will be removed from disk on cleanup.
    GisticQ : double
        A scalar in the range [0, 1] specifying the GISTIC significance
        threshold to use when filtering copy-number events.
    FilterGenes:
        A list of HUGO gene symbols used to further filter significant events.
        Genes that are GISTIC significant but not present on this list will be
        discarded. Can be used with Sanger Cancer Gene Census.
        Default value = None.
    Returns
    -------
    CNVArm : named_tuple
        A named tuple containing the following fields:
        'Symbols' - a numpy array containing chromosome arm symbols of GISTIC
                    significant arm-level copy number events.
        'Description' - comments on the type of alteration either
                        'Amplification', 'Deletion' or
                        'Amplification/Deletion'.
        'CNV' - arm-level log2 ratio copy number values with rows corresponding
                to 'Symbols' and columns corresponding to 'Barcodes'.
        'Barcodes' - numpy array containing TCGA barcodes for samples.
        'Type' - data type - 'CNV-arm' in this case.
        'Release' - release data (string) of GISTIC data from Broad GDAC.
    CNVGene : named_tuple
        A named tuple containing the following fields:
        'Symbols' - a numpy array containing HUGO gene symbols of GISTIC
                    significant gene-level copy number events.
        'Description' - comments on the type of alteration either
                        'Amplification', 'Deletion' or
                        'Amplification/Deletion'.
        'CNV' - gene-level log2 ratio copy number values with rows
                corresponding to 'Symbols' and columns corresponding to
                'Barcodes'.
        'Barcodes' - numpy array containing TCGA barcodes for samples.
        'Type' - data type - 'CNV-arm' in this case.
        'Release' - release data (string) of GISTIC data from Broad GDAC.

    Notes
    -----
    Return values are returned as a namedtuple.
    """

    # make output folder if exists
    if not os.path.isdir(Output):
        os.mkdir(Output)

    # get name of latest firehose data run
    FH = subprocess.Popen([FirehosePath + "firehose_get -r"],
                          stdout=subprocess.PIPE, shell=True)
    (Runs, err) = FH.communicate()
    Latest = [Run for Run in Runs.split("\n") if
              Run.startswith("analyses")][-1]

    # fetch mutation data from firehose and move to output
    FH = subprocess.Popen(["cd " + Output + "; " + FirehosePath +
                           "firehose_get -b -tasks " +
                           "CopyNumber_Gistic2.Level_4 analysis latest " +
                           Disease], stdout=subprocess.PIPE, shell=True)
    (Out, err) = FH.communicate()

    # find *.tar.gz files in download directory
    Files = []
    for root, dirs, files in os.walk(Output + Latest):
        for file in files:
            if file.endswith(".tar.gz"):
                Files.append(root + "/" + file)

    # extract GISTIC files
    MAFZip = [File for File in Files if
              "CopyNumber_Gistic2.Level_4" in File][0]
    Tar = tarfile.open(MAFZip)
    BroadSigFile = [member for member in Tar.getmembers() if
                    "broad_significance_results.txt" in member.name][0]
    BroadValueFile = [member for member in Tar.getmembers() if
                      "broad_values_by_arm.txt" in member.name][0]
    AmpConfFile = [member for member in Tar.getmembers() if
                   "table_amp.conf_99.txt" in member.name][0]
    DelConfFile = [member for member in Tar.getmembers() if
                   "table_del.conf_99.txt" in member.name][0]
    ThreshFile = [member for member in Tar.getmembers() if
                  "all_thresholded.by_genes.txt" in member.name][0]
    BroadSigFile.name = os.path.basename(BroadSigFile.name)
    BroadValueFile.name = os.path.basename(BroadValueFile.name)
    AmpConfFile.name = os.path.basename(AmpConfFile.name)
    DelConfFile.name = os.path.basename(DelConfFile.name)
    ThreshFile.name = os.path.basename(ThreshFile.name)
    Tar.extract(BroadSigFile, path=Output)
    Tar.extract(BroadValueFile, path=Output)
    Tar.extract(AmpConfFile, path=Output)
    Tar.extract(DelConfFile, path=Output)
    Tar.extract(ThreshFile, path=Output)

    # read in arm-level significance results
    TextFile = open(Output + "broad_significance_results.txt", 'r')
    Contents = [line[:-1].split('\t') for line in TextFile]
    ArmIndex = Contents[0].index("Arm")
    AmpIndex = Contents[0].index("Amp q-value")
    DelIndex = Contents[0].index("Del q-value")
    Contents = np.array(Contents[1:])
    Arms = Contents[:, ArmIndex]
    AmpQs = Contents[:, AmpIndex].astype(np.float)
    DelQs = Contents[:, DelIndex].astype(np.float)
    AmpArms = Arms[AmpQs <= GisticQ]
    DelArms = Arms[DelQs <= GisticQ]
    SigArms = list(set(DelArms).union(set(AmpArms)))
    SigArms = [Arm for Arm in SigArms if not((Arm == 'Xp') | (Arm == 'Xq'))]
    SigArms.sort()
    ArmType = []
    for Arm in SigArms:
        if Arm in AmpArms and not (Arm in DelArms):
            ArmType.append('Amplification')
        if not (Arm in AmpArms) and Arm in DelArms:
            ArmType.append('Deletion')
        if Arm in AmpArms and Arm in DelArms:
            ArmType.append('Amplification/Deletion')

    # read in arm-level values and filter out insignificant events
    TextFile = open(Output + "broad_values_by_arm.txt", 'r')
    Contents = np.array([line[:-1].split('\t') for line in TextFile])
    ArmSymbols = Contents[1:, 0]
    ArmBarcodes = Contents[0, 1:]
    ArmCNV = Contents[1:, 1:].astype(np.float)
    Keep = [np.nonzero(ArmSymbols == Arm)[0][0] for Arm in SigArms]
    Keep.sort()
    ArmSymbols = ArmSymbols[Keep]
    ArmCNV = ArmCNV[Keep, :]

    # read in amplification conf file to generate amplification genes
    TextFile = open(Output + "table_amp.conf_99.txt", 'r')
    Contents = [line[:-1].split('\t') for line in TextFile]
    GeneIndex = Contents[0].index("genes_in_region")
    Contents = np.array(Contents[1:])
    GeneLists = Contents[1:, GeneIndex]
    AmpGenes = [List.split(',') for List in GeneLists]
    AmpGenes = [Gene for List in AmpGenes for Gene in List]

    # read in deletion conf file to generate amplification genes
    TextFile = open(Output + "table_del.conf_99.txt", 'r')
    Contents = [line[:-1].split('\t') for line in TextFile]
    GeneIndex = Contents[0].index("genes_in_region")
    Contents = np.array(Contents[1:])
    GeneLists = Contents[1:, GeneIndex]
    DelGenes = [List.split(',') for List in GeneLists]
    DelGenes = [Gene for List in DelGenes for Gene in List]

    # build unified gene list containing amplifications and deletions
    SigGenes = list(set(DelGenes).union(set(AmpGenes)))
    SigGenes = [Gene for Gene in SigGenes if not Gene == '']

    # remove leading and trailing square brackets where present
    SigGenes = [Gene[1:] if Gene[0] == '[' else Gene for Gene in SigGenes]
    SigGenes = [Gene[0:-1] if Gene[-1] == ']' else Gene for Gene in SigGenes]

    # sort gene symbols alphabetically
    SigGenes.sort()

    # determine alteration type for each gene
    GeneType = []
    for Gene in SigGenes:
        if Gene in AmpGenes and not (Gene in DelGenes):
            GeneType.append('Amplification')
        if not (Gene in AmpGenes) and Gene in DelGenes:
            GeneType.append('Deletion')
        if Gene in AmpGenes and Gene in DelGenes:
            GeneType.append('Amplification/Deletion')

    # read in thresholded copy number values and process
    TextFile = open(Output + "all_thresholded.by_genes.txt", 'r')
    Contents = np.array([line[:-1].split('\t') for line in TextFile])
    Symbols = Contents[1:, 0]
    Barcodes = Contents[0, 3:]
    CNV = Contents[1:, 3:].astype(np.float)
    Keep = [np.nonzero(Symbols == Gene)[0][0] for Gene in SigGenes]
    Symbols = Symbols[Keep]
    CNV = CNV[Keep, :]

    # filter out non-protein coding genes and non-HUGO gene symbols
    Symbols, Indices = HUGOFilter(list(Symbols))
    CNV = CNV[Indices, :]
    GeneType = [GeneType[Index] for Index in Indices]

    # filter genes using user provided list
    if FilterGenes is not None:
        Indices = [Symbols.index(Gene) for Gene in Symbols
                   if Gene in FilterGenes]
        Symbols = [Symbols[Index] for Index in Indices]
        CNV = CNV[Indices, :]
        GeneType = [GeneType[Index] for Index in Indices]

    # cleanup files
    os.remove(Output + "broad_significance_results.txt")
    os.remove(Output + "broad_values_by_arm.txt")
    os.remove(Output + "table_amp.conf_99.txt")
    os.remove(Output + "table_del.conf_99.txt")
    os.remove(Output + "all_thresholded.by_genes.txt")
    shutil.rmtree(Output + Latest)

    # build named tuples for outputs
    CNVTuple = namedtuple('CNV', ['Symbols', 'Description', 'CNV', 'Barcodes',
                                  'Type', 'Release'])
    CNVArm = CNVTuple(ArmSymbols, ArmType, ArmCNV, ArmBarcodes, 'CNV-arm',
                      Latest)
    CNVGene = CNVTuple(Symbols, GeneType, CNV, Barcodes, 'CNV-gene',
                       Latest)

    # return outputs
    return CNVArm, CNVGene
