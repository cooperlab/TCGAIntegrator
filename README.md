## TCGAIntegrator
TCGAIntegrator is a tool for building python variables taht contain sample-level integrated views of TCGA projects. Given a disease type, TCGAIntegrator will use Firehose (developed by the Broad Institute) to programmatically access databases at the Broad Institute Genome Data Analysis Center to download and integrate the latest clinical, mutation, arm-level copy number, focal copy number, protein expression and gene/mRNA expression data. User provided thresholds are used with Mutsig2CV and GISTIC2 analyses to filter non-significant mutation and copy number events from the integration. A standard set of Clinical Data Elements (CDEs) are defined and integrated by default to capture basic demographic (age, gender, race) survival (time-to-event and vital status at last followup) and treatment (radiation therapy) information. Additional CDEs and their descriptions can be found here: https://tcga-data.nci.nih.gov/docs/dictionary/. A full list of available projects from Firehose is available here: https://gdac.broadinstitute.org/.

## Dependencies
TCGAIntegrator requires the *requests*, *firebrowse* and *numpy* packages.

## Filtering Copy Number Events
Copy number events identified as significant by GISTIC can number in the thousands. These results can be further filtered against the Sanger Cancer Gene Census (http://cancer.sanger.ac.uk/census/) by providing a path to the Census .csv file that can be downloaded from the Sanger project page.

## Encoding Clinical Variables
The default CDEs define a core set of clinical data that is broadly important in many cancers. Disease-specific analyses should examine the TCGA data and define the relevant CDEs to capture important clinical information for that specific application. TCGAIntegrator will first examine the text values of the clinical data obtained with Firehose and attempt to convert each CDE to floats if possible. CDEs that are categorical in nature will be encoded as a sequence of binary variables with appropriate names generated for each case. For example, a CDE 'gender' with possible values 'male' and 'female' will generate a feature with the symbol 'gender-Is-male' and values '0' (no) or '1' (yes).

##Usage
Build and activate a virtual environment at the command line and install the *firebrowse* and *requests* packages:
```
>virtualenv Integrator
>source ./Integrator/bin/activate
>pip install requests
>pip install firebrowse
```

Import the function "BuildDataset" from the library:

```python
from BuildDataset import BuildDataset
```

Define the output folder where the results will be stored, the disease type, the Mutsig and GISTIC thresholds, and call BuildDatasets to build the archive:

```python
Output = '/home/me/output/'
FirehosePath = None
Disease = 'LGG'
CancerCensusFile = None
MutsigQ = 0.1
GisticQ = 0.1
BuildDataset(Output, FirehosePath, Disease, CancerCensusFile, MutsigQ, GisticQ)
```

Deactivate the virtual environment when the script is finished running:

```
>deactivate
```

## Description of Outputs
The script generates a pickle file containing the following variables:

1. Features - a D x N float numpy array where each column contains the integrated profile of a sample, and each row represents a single clinical/genomic feature.  
2. Symbols - a D-length list of strings containing the feature names. These can be gene symbols, chromosome arms, or CDEs. Each symbol is appended with the feature type ('Clinical', 'Mut', 'CNV', 'CNVArm', 'Protein', or 'mRNA')  
3. SymbolTypes - a D-length list of strings containing stand-alone feature types 'Clinical', 'Mut', 'CNV', 'CNVArm', 'Protein', or 'mRNA'.  
4. Samples - an N-length list of strings containing the TCGA barcodes of each sample TCGA-XX-YYYY-ZZ.  
5. Survival - an N-length float numpy array containing the death or last followup times in days for each sample. These are obtained from the CDEs 'days_to_death', 'days_to_last_followup' and 'vital_status'.  
6. Censored - an N-length float numpy array containing the right-censoring status of each sample. A value of '1' indicates samples where the patient was alive at last followup and a value of '0' indicates uncensored samples where a death even was observed.  

An example of selecting specific features (clinical + protein) can be done using list comprehension and slicing:

```python
Indices = [Index for Index, Type in enumerate(SymbolTypes) if Type in ['Clinical', 'Protein']]
Selected = Features[Indices, :]
```
