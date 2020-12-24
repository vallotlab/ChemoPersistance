# H3K27me3 is a determinant of chemotolerance in triple-negative breast cancer

Scripts to analyze single-cell RNA, single-cell ChIP-seq and bulk ChIP-seq of TNBC
under drug treatment. Follow the steps to reproduce the analysis in the paper '**H3K27me3 is a determinant of chemotolerance in triple-negative breast cancer**' by Marsolier et al.

## 0.0 Setup

In order to re-run the analysis from the paper you must first download this repository.
Then, at the base of the repository, create an "output" and "input" empty folders.
Download the processed data (e.g. count matrices & bigwig files) from GSEXXXXX,
and place it in the "input" folder. In the "input" folder, each kind of data should
be placed in a properly created directory (scRNAseq / scChIPseq / bulk_ChIPseq).
The following hiearchy should be kept:  

```
input/
├── bulk_ChIPseq  
│   ├── MM468  
│   │   ├── BigWigs  
│   │   └── ChIPreChIP  
│   └── PDX  
│       ├── BigWigs  
│       └── Count_Matrices  
├── Packages  
├── scChIPseq  
│   ├── MM468  
│   │   ├── BigWigs  
│   │   ├── Count_Matrices  
│   │   └── Raw_Counts  
│   │       ├── K27  
│   │       └── K4  
│   └── PDX  
│       ├── BigWigs  
│       └── Count_Matrices  
└── scRNAseq  
    ├── MM468  
    │   └── barcoding_output  
    └── PDX  
    
```
Please refer to the scripts if you have doubts where you should place your input 
files.

## 1.0 scRNA
### PDX model
Run the scripts in the order (1.0, 2.0, 3.0).  

### MDA-MB-468 model
There are 3 sub-analyses in this directory: the main analysis ('1.Persister') and
two sub analyses in response to epigenomic drug treatment ('2.UNC' and '3.UNC_5FU').
First run the 1.0 and 2.0 QC scripts.
Then you can run the script of each analyses indepentently by respecting the order
in each given analyses. Note that you must first run the single-cell ChIPseq analyse
before running the script '*5.0_ComparisonChIPseq.R*'.


## 2.0 scChIP


## 3.0 bulk ChIPseq