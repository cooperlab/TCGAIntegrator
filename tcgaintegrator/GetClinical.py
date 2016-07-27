from collections import namedtuple
import numpy as np
import os
import shutil
import subprocess
import tarfile


def GetClinical(Output, FirehosePath, Disease,
                FilterCDEs=['age_at_initial_pathologic_diagnosis',
                            'gender', 'histological_type', 'pathologic_stage',
                            'pathologic_m', 'pathologic_n',
                            'pathologic_t', 'race', 'radiation_therapy']):
    """Generates variables containing clinical data values describing patient
    demographics, disease phenotypes and treatment in addition to survival and
    vital status at last followup. Uses Firebrowse, a tool from the Broad
    Genome Data Analysis Center to download clinical data values from the Broad
    Institute servers. Automatically cleans up results on completion.

    Parameters
    ----------
    Output : string
        Path to be used for temporary downloading and unzipping clinical data
        files. Downloads and extracted files will be removed from disk on
        cleanup.
    FirehosePath : string
        Path to firehose_get executable.
    Disease : string
        Dataset code to generate protein expression profiles for. Can be
        obtained using firehose_get -c.
    FilterCDEs : list
        List of strings defining the clinical data elements to return. Survival
        defining terms including 'days_to_death', 'days_to_last_followup', and
        'vital_status' are obtained by default and do not need to be defined
        here. Default CDEs are selected as those defined for a broad set of
        diseases and clinically-relevant.
        Default value = ['age_at_initial_pathologic_diagnosis',
                         'gender', 'histological_type', 'pathologic_stage',
                         'pathologic_m', 'pathologic_n',
                         'pathologic_t', 'race', 'radiation_therapy']

    Returns
    -------
    Clinical : named_tuple
        A named tuple containing the following fields:
        'CDEs' - a numpy array containing clinical data element field names.
        'Values' - numpy array of clinical data values as strings where rows
                       correspond to 'Symbols' and columns corresponding to
                       'Barcodes'.
        'Class' - a numpy array containing strings describing the types of
                  data in the rows of 'Values'. Possible values include
                  'numeric', 'ordinal', or 'binary'.
        'Barcodes' - numpy array containing TCGA barcodes for samples.
        'Type' - data type - 'Clinical' in this case.
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

    # fetch clinical data from firehose and move to output
    FH = subprocess.Popen(["cd " + Output + "; " + FirehosePath +
                           "firehose_get -b -tasks " +
                           "Clinical_Pick_Tier1.Level_4 data latest " +
                           Disease], stdout=subprocess.PIPE, shell=True)
    (Out, err) = FH.communicate()

    # find *.tar.gz files in download directory
    Files = []
    for root, dirs, files in os.walk(Output + Latest):
        for file in files:
            if file.endswith(".tar.gz"):
                Files.append(root + "/" + file)

    # extract clinical data "All_CDEs.txt" file
    ClinicalZip = [File for File in Files if
                   "Clinical_Pick_Tier1.Level_4" in File]
    Tar = tarfile.open(ClinicalZip[0])
    ClinicalFile = [member for member in Tar.getmembers() if
                    "All_CDEs.txt" in member.name]
    ClinicalFile[0].name = os.path.basename(ClinicalFile[0].name)
    Tar.extract(ClinicalFile[0], path=Output)
    Tar.close()

    # extract clinical symbols, values in string form
    TextFile = open(Output + ClinicalFile[0].name, 'r')
    Contents = np.array([line[:-1].split('\t') for line in TextFile])
    TextFile.close()
    Barcodes = list(Contents[0, 1:])
    CDEs = list(Contents[2:, 0])
    Values = Contents[2:, 1:]

    # convert barcodes to uppercase characters
    Barcodes = [Barcode.upper() for Barcode in Barcodes]

    # replace missing values 'NA' with conversion-compatible 'NaN'
    Values[Values == 'NA'] = 'NaN'

    # strip out required survival data
    D2D = Values[CDEs.index('days_to_death'), :].astype(np.float)
    D2LF = Values[CDEs.index('days_to_last_followup'), :].astype(np.float)
    VS = Values[CDEs.index('vital_status'), :]

    # convert survival data to numeric values
    Censored = (VS == 'alive').astype(np.float)
    Censored[VS == 'NaN'] = np.NaN
    Survival = np.NaN * np.zeros(Censored.shape)
    Survival[Censored == 1] = D2LF[Censored == 1]
    Survival[Censored == 0] = D2D[Censored == 0]

    # remove survival related fields from output features
    Values = np.delete(Values, CDEs.index('days_to_death'), axis=0)
    CDEs = [CDE for CDE in CDEs if CDE != 'days_to_death']
    Values = np.delete(Values, CDEs.index('days_to_last_followup'), axis=0)
    CDEs = [CDE for CDE in CDEs if CDE != 'days_to_last_followup']
    Values = np.delete(Values, CDEs.index('vital_status'), axis=0)
    CDEs = [CDE for CDE in CDEs if CDE != 'vital_status']

    # strip out user-defined clinical fields for features
    Indices = [CDEs.index(CDE) for CDE in FilterCDEs if CDE in CDEs]
    CDEs = [CDE for CDE in FilterCDEs if CDE in CDEs]
    Values = Values[Indices, :]

    # iterate over CDEs, converting to numeric values
    Encoded = []
    Names = []
    for i in range(len(CDEs)):
        Numeric = [_CheckNumeric(Value) for Value in Values[i, :]]
        if(sum(Numeric) == len(Numeric)):  # feature is numeric - convert float
            Encoded.append([float(Value) for Value in Values[i, :]])
            Names.append(CDEs[i])  # feature is numeric - convert to float
        else:  # feature is categorical, convert to a series of binary values
            Dict = list(set(list(Values[i, Values[i, :] != 'NaN'])))
            if(len(Dict) == 2):  # binary feature - represent with 1 feature
                Encode = [1 if Value == Dict[0] else 0
                          for Value in Values[i, :]]
                for k in range(len(Encode)):
                    if(Values[i, k] == 'NaN'):
                        Encode[k] = np.NaN
                Name = CDEs[i] + '-Is-' + Dict[0]
                Encoded.append(Encode)
                Names.append(Name)
            else:
                for j in range(len(Dict)):
                    Encode = [1 if Value == Dict[j] else 0
                              for Value in Values[i, :]]
                    for k in range(len(Encode)):
                        if(Values[i, k] == 'NaN'):
                            Encode[k] = np.NaN
                    Name = CDEs[i] + '-Is-' + Dict[j]
                    Encoded.append(Encode)
                    Names.append(Name)

    # convert clinical values to floating-point numpy array
    Values = np.array(Encoded)

    # cleanup files
    os.remove(Output + ClinicalFile[0].name)
    shutil.rmtree(Output + Latest)

    # build named tuples for outputs
    ClinicalTuple = namedtuple('Clinical', ['Survival', 'Censored', 'CDEs',
                                            'Values', 'Barcodes', 'Type',
                                            'Release'])
    Clinical = ClinicalTuple(Survival, Censored, Names, Values, Barcodes,
                             'Clinical', Latest)

    # return outputs
    return Clinical


def _CheckNumeric(s):
    try:
        float(s)
        return 1
    except ValueError:
        return 0
