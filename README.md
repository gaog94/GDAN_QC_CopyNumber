# GDAN_QC_CopyNumber
Perform Copy Number Comparison between TCGA Legacy and Harmonized Data (for Genomic Data Analysis Network aka GDAN)

Description
-----------
As part of the Genomic Data Analysis Network (GDAN), the Quality Control (QC) working group was tasked with performing a comprehensive analysis of the TCGA hg38 "Harmonized" data housed at the NCI's Genomic Data Commons (GDC) to the "Legacy" hg19 data upon which the majority of the TCGA marker papers were based. This section looks specifically at copy number Affy SNP6.0 data, both at the individual gene level and at the level of focal peaks and driver events.

In brief, we use GISTIC2.0 to analyze tumor copy number profiles from each TCGA disease type, comparing the output from hg19-aligned profiles against those from hg38-aligned profiles. In each scenario, germline-subtracted copy number profiles were generated from Affymetrix Genome-Wide Human SNP6.0 arrays run on tumor samples and corresponding blood normal samples (where available). Pre-computed GISTIC analyses of the hg19-aligned data are available at [firebrowse](https://www.firebrowse.org/), while both the legacy hg19-aligned and the updated hg38-aligned copy number segmentation profiles are available for download from [the GDC](https://portal.gdc.cancer.gov/) and from [FireCloud](https://portal.firecloud.org).

We subsequently analyzed the outputs of these GISTIC2.0 runs (compiled in separate directories named by cancer type) using the scripts in this repository. Here we detail instructions for running these scripts.

To run the raw gene-level values comparison:
```
python scripts/Raw_Call_Comparison.py /path/to/hg38/gistic/results/ /path/to/hg19/gistic/results/
```

To run the thresholded gene-level values comparison:
```
python scripts/Thresholded_Call_Comparison.py /path/to/hg38/gistic/results/ /path/to/hg19/gistic/results/ /path/to/aneuploidy/focal/alteration/file.txt
```

To visualize overlaps in genes called in significant focal peaks:
```
python scripts/Focal_Peak_Gene_Overlap.py /path/to/hg38/gistic/results/ /path/to/hg19/gistic/results/
```

To print a quick summary of which documented copy number driver genes (this list was manually compiled from the TCGA marker papers) are still conserved/called in legacy and harmonized runs:
```
python scripts/Driver_Conservation.py ../Documented_Driver_Differences.txt
```
