from collections import namedtuple
from HUGOFilter import HUGOFilter
import numpy as np
import os
import shutil
import subprocess
import tarfile


def GetGeneExpression(FirehosePath, Disease, Output):
    """Generates variables containing gene expression values from the
    platform. Uses Firebrowse, a tool from the Broad Genome Data
    Analysis Center to download gene expression array values from the Broad
    Institute servers. Automatically cleans up results on completion.

    Parameters
    ----------
    FirehosePath : string
        Path to firehose_get executable.
    Disease : string
        Dataset code to generate protein expression profiles for. Can be
        obtained using firehose_get -c.
    Output : string
        Path to be used for temporary downloading and unzipping mRNA
        files. Downloads and extracted files will be removed from disk on
        cleanup.

    Returns
    -------
    mRNA : named_tuple
        A named tuple containing the following fields:
        'Symbols' - a numpy array containing gene symbols from protein
                    expression arrays. Each entry may contain multiple symbols
                    as this assay uses antibodies to measure expression.
        'Description' - Ensemble IDs associated with 'Symbols'.
        'Expression' - numpy array of gene expression values where rows
                       correspond to 'Symbols' and columns corresponding to
                       'Barcodes'.
        'Barcodes' - numpy array containing TCGA barcodes for samples.
        'Type' - data type - 'GeneExpression' in this case.
        'Release' - release data (string) of data from Broad GDAC.

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
                           "Level_3__RSEM_genes_normalized__data.Level_3 " +
                           "stddata latest " + Disease],
                          stdout=subprocess.PIPE, shell=True)
    (Out, err) = FH.communicate()

    # find *.tar.gz files in download directory
    Files = []
    for root, dirs, files in os.walk(Output + Latest):
        for file in files:
            if file.endswith(".tar.gz"):
                Files.append(root + "/" + file)

    # extract gene expression values file
    GEZip = [File for File in Files if
             "Level_3__RSEM_genes_normalized__data.Level_3" in File]
    Tar = tarfile.open(GEZip[0])
    GEFile = [member for member in Tar.getmembers() if member.name.find(
              "Level_3__RSEM_genes_normalized__data.data.txt") != -1]
    GEFile[0].name = os.path.basename(GEFile[0].name)
    Tar.extract(GEFile[0], path=Output)
    Tar.close()

    # extract gene symbols, antibodies and gene expression values
    TextFile = open(Output + GEFile[0].name, 'r')
    Contents = np.array([line[:-1].split('\t') for line in TextFile])
    TextFile.close()
    Barcodes = list(Contents[0, 1:])
    Symbols = list(Contents[2:, 0])
    Values = Contents[2:, 1:]
    Values[Values == "NA"] = "nan"
    Expression = np.zeros(Values.shape)
    for i in range(Values.shape[0]):
        Expression[i, :] = np.genfromtxt(Values[i, :])
    Description = [Symbol.split('|')[1].strip() for Symbol in Symbols]
    Symbols = [Symbol.split('|')[0].strip() for Symbol in Symbols]

    # filter out non-protein coding genes and non-HUGO gene symbols
    Symbols, Indices = HUGOFilter(list(Symbols))
    Expression = Expression[Indices, :]
    Description = [Description[Index] for Index in Indices]

    # cleanup files
    os.remove(Output + GEFile[0].name)
    shutil.rmtree(Output + Latest)

    # build named tuples for outputs
    Tuple = namedtuple('Gene', ['Symbols', 'Description',
                                'Expression', 'Barcodes',
                                'Type', 'Release'])
    mRNA = Tuple(Symbols, Description, Expression, Barcodes,
                 'GeneExpression', Latest)

    # return outputs
    return mRNA
