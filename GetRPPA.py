from collections import namedtuple
import numpy as np
import os
import shutil
import subprocess
import tarfile


def GetRPPA(FirehosePath, Disease, Output):
    """Generates variables containing protein expression values from the RPPA
    platform. Uses Firebrowse, a tool from the Broad Genome Data
    Analysis Center to download RPPA array values from the Broad Institute
    servers. Automatically cleans up results on completion.

    Parameters
    ----------
    FirehosePath : string
        Path to firehose_get executable.
    Disease : string
        Dataset code to generate protein expression profiles for. Can be
        obtained using firehose_get -c.
    Output : string
        Path to be used for temporary downloading and unzipping RPPA
        files. Downloads and extracted files will be removed from disk on
        cleanup.

    Returns
    -------
    Protein : named_tuple
        A named tuple containing the following fields:
        'Symbols' - a numpy array containing gene symbols from protein
                    expression arrays. Each entry may contain multiple symbols
                    as this assay uses antibodies to measure expression.
        'Description' - comments on the antibody associated with 'Symbols'.
        'Expression' - numpy array of protein expression values where rows
                       correspond to 'Symbols' and columns corresponding to
                       'Barcodes'.
        'Barcodes' - numpy array containing TCGA barcodes for samples.
        'Type' - data type - 'ProteinExpression' in this case.
        'Release' - release data (string) of RPPA data from Broad GDAC.

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
    Latest = [Run for Run in Runs.split("\n") if Run.startswith("stddata")][-1]

    # fetch RPPA data from firehose and move to output
    FH = subprocess.Popen(["cd " + Output + "; " + FirehosePath +
                           "firehose_get -b -tasks " +
                           "RPPA_AnnotateWithGene.Level_3 stddata latest " +
                           Disease], stdout=subprocess.PIPE, shell=True)
    (Out, err) = FH.communicate()

    # find *.tar.gz files in download directory
    Files = []
    for root, dirs, files in os.walk(Output + Latest):
        for file in files:
            if file.endswith(".tar.gz"):
                Files.append(root + "/" + file)

    # extract Disease.rppa.txt file
    RPPAZip = [File for File in Files if
               "RPPA_AnnotateWithGene.Level_3" in File]
    Tar = tarfile.open(RPPAZip[0])
    RPPAFile = [member for member in Tar.getmembers() if Disease +
                ".rppa.txt" in member.name]
    RPPAFile[0].name = os.path.basename(RPPAFile[0].name)
    Tar.extract(RPPAFile[0], path=Output)
    Tar.close()

    # extract gene symbols, antibodies and RPPA expression values
    TextFile = open(Output + RPPAFile[0].name, 'r')
    Contents = np.array([line[:-1].split('\t') for line in TextFile])
    TextFile.close()
    Barcodes = list(Contents[0, 1:])
    Symbols = list(Contents[1:, 0])
    Values = Contents[1:, 1:]
    Values[Values == "NA"] = "nan"
    Expression = np.zeros(Values.shape)
    for i in range(Values.shape[0]):
        Expression[i, :] = np.genfromtxt(Values[i, :])
    Description = [Symbol.split('|')[1].strip() for Symbol in Symbols]
    Symbols = [Symbol.split('|')[0].strip() for Symbol in Symbols]

    # cleanup files
    os.remove(Output + RPPAFile[0].name)
    shutil.rmtree(Output + Latest)

    # build named tuples for outputs
    ProteinTuple = namedtuple('Protein', ['Symbols', 'Description',
                                          'Expression', 'Barcodes',
                                          'Type', 'Release'])
    Protein = ProteinTuple(Symbols, Description, Expression, Barcodes,
                           'ProteinExpression', Latest)

    # return outputs
    return Protein
