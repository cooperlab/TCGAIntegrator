from firebrowse import fbget
from .GetClinical import GetClinical
from .GetCopyNumber import GetCopyNumber
from .GetGeneExpression import GetGeneExpression
from .GetMutations import GetMutations
from .GetRPPA import GetRPPA
import numpy as np
import os
import scipy.io as io
import subprocess
import sys
import time
import zipfile


def BuildDataset(Output, FirehosePath=None, Disease=None,
                 CancerCensusFile=None, MutsigQ=0.1, Raw=False, GisticQ=0.25,
                 FilterCDEs=['age_at_initial_pathologic_diagnosis',
                             'days_to_death', 'days_to_last_followup',
                             'gender', 'histological_type', 'pathologic_stage',
                             'pathologic_m', 'pathologic_n',
                             'pathologic_t', 'race', 'radiation_therapy',
                             'vital_status'], SampleCodes=[1, 2]):
    """Generates TCGA data in python formats for mRNA expression, protein
    expression, copy number, mutation and clinical platforms. All data is
    automatically downloaded and curated from the Broad Institute GDAC using
    their firehose_get tool. Uses results from Mutsig2CV and GISTIC algorithms
    to filter mutation and copy number data to cancer-relevant genetic events.

    Parameters
    ----------
    Output : string
        Path to be used for generating outputs. Temporary downloading and
        unzipping of files will also happen here. Downloads and extracted files
        will be removed from disk on cleanup.
    FirehosePath : string
        Path to firehose_get executable. If not provided the executable will
        be downloaded to the folder 'Output' using the command:
        'wget http://gdac.broadinstitute.org/runs/code/firehose_get_latest.zip
        -P Output'
        and unpacked in that location.
    Disease : string
        Dataset code to generate protein expression profiles for. Can be
        obtained using firehose_get -c. If not provided, all datasets will be
        built. Default value = None.
    CancerCensusFile : string
        Path and filename for the Sanger Cancer Gene Census .csv file, obtained
        from http://cancer.sanger.ac.uk/census/. Used for filtering copy
        number events identified as significant by GISTIC.
        Default value = None.
    MutsigQ : double
        A scalar in the range [0, 1] specifying the Mutsig2CV significance
        threshold to use when filtering somatic mutation events.
    Raw : bool
        Flag indicating whether to use raw mutation calls, or packaged mutation
        calls. Packaged calls may have fewer samples. Raw calls may have more
        samples but are unavailable for some projects.        
    GisticQ : double
        A scalar in the range [0, 1] specifying the GISTIC significance
        threshold to use when filtering copy-number events.
    FilterCDEs : list
        List of strings defining the clinical data elements to return. Default
        CDEs are selected as those defined for a broad set of diseases and
        clinically-relevant.
        Default value = ['age_at_initial_pathologic_diagnosis',
                         'gender', 'histological_type', 'pathologic_stage',
                         'pathologic_m', 'pathologic_n',
                         'pathologic_t', 'race', 'radiation_therapy']
    SamplesCodes : list
        List of integer codes identifying sample types to keep. The full list
        of sample type codes can be found at https://tcga-data.nci.nih.gov.
        The default value corresponds to 1 - primary tumor and 2 - recurrent
        tumor. Default value = [1, 2].

    Returns
    -------
    Generates a sequence of files 'Protein', 'mRNA', 'Mutation', 'CNV' and
    'Clinical' for the disease of interest in the folder 'Output'. If 'Disease'
    was selected as 'None' and data for all diseases are generated, then
    these data will be organized into subdirectories under 'Output' by disease.
    """

    # add trailing slash to 'Output' if present
    if Output[-1] != '/':
        Output = Output + "/"

    # download firehose_get if path not provided
    if FirehosePath is None:
        
        # set flag to remove executable later
        Remove = True

        # make output directory of not available
        if not os.path.isdir(Output):
            os.mkdir(Output)

        # wget firehose_get binary to Output folder, unzip and delete .zip
        FH = subprocess.Popen("wget http://gdac.broadinstitute.org/runs/code/"
                              "firehose_get_latest.zip -P " + Output,
                              stdout=subprocess.PIPE, shell=True)
        (Runs, err) = FH.communicate()
        Zip = zipfile.ZipFile(Output + "firehose_get_latest.zip")
        Zip.extractall(Output)
        os.remove(Output + "firehose_get_latest.zip")
        os.chmod(Output + "firehose_get", 0755)

        # set 'FirehosePath' to Output
        FirehosePath = Output

    # check if disease was provided
    if Disease is None:

        # get list of available diseases
        Chars = fbget.cohorts()
        Cohorts = [row.split("\t") for row in Chars.split("\n")]
        Diseases = [str(Cohorts[i][0]) for i in range(1, len(Cohorts)-1)]

        # generate prefix necessary for creating subdirectories
        Prefixes = [Cohort + "/" for Cohort in Diseases]

    else:

        # process a single disease
        Diseases = [Disease]

        # output will be put directly into 'Output'
        Prefixes = ['']

    # extract Sanger Census genes if path is provided
    if CancerCensusFile is not None:
        File = open(CancerCensusFile, 'r')
        CancerGenes = [line[:-1].split(',')[0] for line in File]
    else:
        CancerGenes = None

    # iterate over each disease type, generating output and updating console
    for CohortIndex, Cohort in enumerate(Diseases):

        # write disease type to console
        sys.stdout.write("Processing " + Cohort + ", Disease " +
                         str(CohortIndex+1) + " of " + str(len(Diseases)) +
                         "\n")
        sys.stdout.write("\tOutput will be generated in " + Output + "\n")

        # generate clinical data
        sys.stdout.write("\tClinical - generating data...")
        Start = time.time()
        Clinical = GetClinical(Output + Prefixes[CohortIndex], FirehosePath,
                               Cohort, FilterCDEs)
        sys.stdout.write(" done in " + str(time.time()-Start) +
                         " seconds.\n")

        # generate mutation
        sys.stdout.write("\tMutations - generating data...")
        Start = time.time()
        Mutations = GetMutations(Output + Prefixes[CohortIndex], FirehosePath,
                                 Cohort, MutsigQ, Raw)
        sys.stdout.write(" done in " + str(time.time()-Start) +
                         " seconds.\n")

        # generate copy number
        sys.stdout.write("\tCopy Number - generating data...")
        Start = time.time()
        (CNVArm, CNVGene) = GetCopyNumber(Output + Prefixes[CohortIndex],
                                          FirehosePath, Cohort, GisticQ,
                                          CancerGenes)
        sys.stdout.write(" done in " + str(time.time()-Start) +
                         " seconds.\n")

        # generate RPPA
        sys.stdout.write("\tProtein Expression - generating data,")
        Start = time.time()
        Protein = GetRPPA(Output + Prefixes[CohortIndex], FirehosePath, Cohort)
        sys.stdout.write(" done in " + str(time.time()-Start) +
                         " seconds.\n")

        # generate gene expression
        sys.stdout.write("\tmRNA expression - generating data,")
        Start = time.time()
        mRNA = GetGeneExpression(Output + Prefixes[CohortIndex], FirehosePath,
                                 Cohort)
        sys.stdout.write(" done in " + str(time.time()-Start) +
                         " seconds.\n")

        # build feature names list
        ClinicalSymbols = [CDE + "_Clinical" for CDE in Clinical.CDEs]
        MutationSymbols = [Symbol + "_Mut" for Symbol in Mutations.Symbols]
        CNVGeneSymbols = [Symbol + "_CNV" for Symbol in CNVGene.Symbols]
        CNVArmSymbols = [Symbol + "_CNVArm" for Symbol in CNVArm.Symbols]
        ProteinSymbols = [Symbol + "_Protein" for Symbol in Protein.Symbols]
        mRNASymbols = [Symbol + "_mRNA" for Symbol in mRNA.Symbols]
        Symbols = ClinicalSymbols + MutationSymbols + CNVGeneSymbols +\
            CNVArmSymbols + ProteinSymbols + mRNASymbols

        # build feature types list
        ClinicalTypes = ["Clinical" for CDE in Clinical.CDEs]
        MutationTypes = ["Mutation" for Symbol in Mutations.Symbols]
        CNVGeneTypes = ["CNVGene" for Symbol in CNVGene.Symbols]
        CNVArmTypes = ["CNVArm" for Symbol in CNVArm.Symbols]
        ProteinTypes = ["Protein" for Symbol in Protein.Symbols]
        mRNATypes = ["mRNA" for Symbol in mRNA.Symbols]
        SymbolTypes = ClinicalTypes + MutationTypes + CNVGeneTypes +\
            CNVArmTypes + ProteinTypes + mRNATypes

        # build comprehensive sample list for molecular types, get tissue code
        MutationSamples = [Barcode[0:15] for Barcode in Mutations.Barcodes]
        CNVGeneSamples = [Barcode[0:15] for Barcode in CNVGene.Barcodes]
        CNVArmSamples = [Barcode[0:15] for Barcode in CNVArm.Barcodes]
        ProteinSamples = [Barcode[0:15] for Barcode in Protein.Barcodes]
        mRNASamples = [Barcode[0:15] for Barcode in mRNA.Barcodes]
        Samples = list(set(MutationSamples + CNVGeneSamples +
                           CNVArmSamples + ProteinSamples + mRNASamples))
        Samples.sort()
        SampleTypes = [int(Sample[13:15]) for Sample in Samples]

        # filter out non-primary tumor tissue types
        Samples = [Sample for Index, Sample in enumerate(Samples)
                   if SampleTypes[Index] in SampleCodes]

        # map clinical samples with short barcodes to 'Samples'
        ClinicalSamples = Clinical.Barcodes
        ClinicalMapped = np.NaN * np.ones((len(ClinicalSymbols), len(Samples)))
        Survival = np.NaN * np.ones((len(Samples)))
        Censored = np.NaN * np.ones((len(Samples)))
        for Current, ClinicalSample in enumerate(ClinicalSamples):
            Indices = [Ind for Ind, Sample in enumerate(Samples)
                       if Sample[0:12] == ClinicalSample]
            ClinicalMapped[:, Indices] = \
                Clinical.Values[:, Current, np.newaxis]
            Survival[Indices] = Clinical.Survival[Current]
            Censored[Indices] = Clinical.Censored[Current]
        AvailableClinical = ['Yes' if Sample[0:12] in ClinicalSamples else 'No'
                             for Sample in Samples]
        AvailableClinical = np.array(AvailableClinical, dtype=object)

        # reshape arrays from mutation data to match order, size of 'Samples'
        Indices = [Samples.index(Sample) for Sample in MutationSamples if
                   Sample in Samples]
        Mapped = [Index for Index, Sample in enumerate(MutationSamples) if
                  Sample in Samples]
        MutationsMapped = np.NaN * np.ones((len(MutationSymbols),
                                            len(Samples)))
        MutationsMapped[:, Indices] = Mutations.Binary[:, Mapped]
        AvailableMutation = ['Yes' if Sample in MutationSamples else 'No'
                             for Sample in Samples]
        AvailableMutation = np.array(AvailableMutation, dtype=object)

        # reshape arrays from CNV data to match order, size of 'Samples'
        Indices = [Samples.index(Sample) for Sample in CNVGeneSamples if
                   Sample in Samples]
        Mapped = [Index for Index, Sample in enumerate(CNVGeneSamples) if
                  Sample in Samples]
        CNVGeneMapped = np.NaN * np.ones((len(CNVGeneSymbols), len(Samples)))
        CNVGeneMapped[:, Indices] = CNVGene.CNV[:, Mapped]
        AvailableCNV = ['Yes' if Sample in CNVGeneSamples else 'No'
                        for Sample in Samples]
        AvailableCNV = np.array(AvailableCNV, dtype=object)
        Indices = [Samples.index(Sample) for Sample in CNVArmSamples if
                   Sample in Samples]
        Mapped = [Index for Index, Sample in enumerate(CNVArmSamples) if
                  Sample in Samples]
        CNVArmMapped = np.NaN * np.ones((len(CNVArmSymbols), len(Samples)))
        CNVArmMapped[:, Indices] = CNVArm.CNV[:, Mapped]

        # reshape arrays from Protein data to match order, size of 'Samples'
        Indices = [Samples.index(Sample) for Sample in ProteinSamples if
                   Sample in Samples]
        Mapped = [Index for Index, Sample in enumerate(ProteinSamples) if
                  Sample in Samples]
        ProteinMapped = np.NaN * np.ones((len(ProteinSymbols), len(Samples)))
        ProteinMapped[:, Indices] = Protein.Expression[:, Mapped]
        AvailableProtein = ['Yes' if Sample in ProteinSamples else 'No'
                            for Sample in Samples]
        AvailableProtein = np.array(AvailableProtein, dtype=object)
                            
        # reshape arrays from Protein data to match order, size of 'Samples'
        Indices = [Samples.index(Sample) for Sample in mRNASamples if
                   Sample in Samples]
        Mapped = [Index for Index, Sample in enumerate(mRNASamples) if
                  Sample in Samples]
        mRNAMapped = np.NaN * np.ones((len(mRNASymbols), len(Samples)))
        mRNAMapped[:, Indices] = mRNA.Expression[:, Mapped]
        AvailablemRNA = ['Yes' if Sample in mRNASamples else 'No'
                        for Sample in Samples]
        AvailablemRNA = np.array(AvailablemRNA, dtype=object)
        
        # stack into master table
        Features = np.vstack((ClinicalMapped, MutationsMapped, CNVGeneMapped,
                              CNVArmMapped, ProteinMapped, mRNAMapped))

        # write outputs to disk
        io.savemat(Output + Prefixes[CohortIndex] + Cohort + ".Data.mat",
                   {'Symbols': Symbols, 'SymbolTypes':SymbolTypes,
                    'Samples':Samples, 'Features':Features,
                    'Survival':Survival, 'Censored':Censored,
                    'AvailableClinical':AvailableClinical,
                    'AvailableCNV':AvailableCNV,
                    'AvailableMutation':AvailableMutation,
                    'AvailableProtein':AvailableProtein,
                    'AvailablemRNA':AvailablemRNA})
        
        # cleanup firehose_get executable if not provided
        if Remove:
            os.remove(Output + "firehose_get")

    return
