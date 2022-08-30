# IDIA

A IDIA tool to generate pseudo-spectra from data independent acquisition (DIA) data mass spectrometry-based proteomics data. The critical innovation of IDAT is that IDIA integrates two isotopic trace detection strategies and employs B-spline and Gaussian filters to help extract high-quality pseudo-spectra from the complex DIA data.

## Setup

### Dependency

- java version 11
- OpenJDK 11.0.6

### Requirement

- 16 GB available RAM or above(based on your data)

## User Manual

Download the source code and import them into Eclipse. Export executable Java JAR file. Run JAR file as:

```
java -jar -Xmx16G IDIA.jar mzMXL_file idia.params
```

mzXML_file: mzXML format mass spectrum data

idia.params: the parameters of IDIA

## Toy example of IDIA

### Data Source of the Toy Example

Download the data via the PRIDE repository PXD001587. The UPS2 DIA data files:

18186_REP2_4pmol_UPS2_SWATH_1.wiff

18186_REP2_4pmol_UPS2_SWATH_1.wiff.scan

### Converting Raw Data

#### First step

Download the AB MS Data Converter: <https://www.sciex.com/form-pages/sw-downloads-form?d=sciex_ms_data_converter_V1.3.1.zip&asset=software&softwareProduct=MS%20Data%20Converter%201.3.1>

The manual of AB MS Data Converter: <https://download.sciex.com/sciex-ms-data-converter-user-guide-en.pdf?_ga=2.238304643.2059936003.1661275969-1581161357.1661016015>
Command line syntax:

```
AB_SCIEX_MS_Converter WIFF "18186_REP2_4pmol_UPS2_SWATH_1.wiff" -profile MZML "18186_REP2_4pmol_UPS2_SWATH_1.mzML"
```

#### Second step

You convert files from mzML format to mzXML formate by msconvert.exe from ProteoWizard using default parameters.  The link of ProteoWizard: <https://proteowizard.sourceforge.io/download.html>

### Run IDIA

Export executable Java JAR file (IDIA.jar) from Eclipse.

Run IDIA as the command line:

```
java -jar -Xmx16G IDIA.jar 18186_REP2_4pmol_UPS2_SWATH_1.mzXML idia.se_params
```

The results:

18186_REP2_4pmol_UPS2_SWATH_1_Q1.mgf and 18186_REP2_4pmol_UPS2_SWATH_1_Q2.mgf
18186_REP2_4pmol_UPS2_SWATH_1_Q1.mzML and 18186_REP2_4pmol_UPS2_SWATH_1_Q2.mzML

There are two types of results (mfg and mzML format). Q1 and Q2 are the different quality of the precursor isotopic envelopes.

### Peptide and Protein Identification

Using all generated pseudo-spectra files with mzML format to identified peptides and proteins by some popular tools such as COMET, PeptideProphet, iProphet, ProteinProphet and so on.

## Feedback

If you have any feedback, please reach out to us at Xuan.Guo@unt.edu and jianchengli@my.unt.edu.
