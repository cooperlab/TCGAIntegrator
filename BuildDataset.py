from firebrowse import fbget
import os
import pickle
import subprocess
import sys
import timeit
import zipfile


def BuildDataset(Output, FirehosePath=None, Disease=None, 
                 MutsigQ=0.1, GisticQ=0.25):
    """Generates TCGA data in python formats for mRNA expression, protein
    expression, copy number, mutation and clinical platforms. All data is
    automatically downloaded and curated from the Broad Institute GDAC using
    their firehose_get tool. Uses results from Mutsig2CV and GISTIC algorithms
    to filter mutation and copy number data to cancer-relevant genetic events.

    Parameters
    ----------
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
    MutsigQ : double
        A scalar in the range [0, 1] specifying the Mutsig2CV significance
        threshold to use when filtering somatic mutation events.
    GisticQ : double
        A scalar in the range [0, 1] specifying the GISTIC significance
        threshold to use when filtering copy-number events.
    Output : string
        Path to be used for generating outputs. Temporary downloading and
        unzipping of files will also happen here. Downloads and extracted files
        will be removed from disk on cleanup.

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

        # make output directory of not available
        if not os.path.isdir(Output):
            os.mkdir(Output)

        # wget firehose_get binary to Output folder, unzip and delete .zip
        Location = "http://gdac.broadinstitute.org/runs/code/firehose_get_"
        "latest.zip"
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
        DiseaseDescription = [str(Cohorts[i][1]) for i in
                              range(1, len(Cohorts)-1)]

        # generate prefix necessary for creating subdirectories
        Prefixes = [Cohort + "/" for Cohort in Diseases]

    else:

        # process a single disease
        Diseases = [Disease]

        # output will be put directly into 'Output'
        Prefixes = ['']

    # iterate over each disease type, generating output and updating console
    for Index, Cohort in enumerate(Diseases):

        # write disease type to console
        sys.stdout.write("Processing " + Cohort + ", Disease " + str(Index+1) +
                         " of " + str(len(Diseases)) + "\n")
        sys.stdout.write("\tOutput will be generated in " + Output + "\n")

        # generate mutation
        sys.stdout.write("\tMutations - generating data...")
        Start = timeit.timeit()
        Mutations = GetMutations(FirehosePath, MutsigQ, Cohort,
                                 Output + Prefixes[Index])
        sys.stdout.write(" done in " + str(timeit.timeit()-Start) +
                         " seconds.\n")

        # generate copy number
        sys.stdout.write("\tCopy Number - generating data...")
        Start = timeit.timeit()
        (CNVArm, CNVGene) = GetCopyNumber(FirehosePath, GisticQ, Cohort,
                                          Output + Prefixes[Index])
        sys.stdout.write(" done in " + str(timeit.timeit()-Start) +
                         " seconds.\n")

        # generate RPPA
        sys.stdout.write("\tProtein Expression - generating data,")
        Start = timeit.timeit()
        Protein = GetRPPA(FirehosePath, Cohort, Output + Prefixes[Index])
        sys.stdout.write(" done in " + str(timeit.timeit()-Start) +
                         " seconds.\n")

        # genearte gene expression
        sys.stdout.write("\tmRNA expression - generating data,")
        Start = timeit.timeit()
        mRNA = GetGeneExpression(FirehosePath, Cohort,
                                 Output + Prefixes[Index])
        sys.stdout.write(" done in " + str(timeit.timeit()-Start) +
                         " seconds.\n")

        # merge outputs and write outputs to disk
        

    return
