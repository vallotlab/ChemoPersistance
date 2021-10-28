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
│   ├── Human
│   │   ├── BigWigs
│   │   └── Count_Matrices
│   ├── InVitro
│   │   ├── bulk_ChIPseq
│   │   │   ├── BigWigs
│   │   │   └── Count_Matrices
│   │   └── ChIPreChIP
│   │       ├── BigWigs_Compare
│   │       │   └── Raw
│   │       ├── BigWigs_Coverage
│   │       └── Peaks
│   └── PDX
│       └── ChIPreChIP
│           ├── BigWigs_Compare
│           │   └── Raw
│           ├── BigWigs_Coverage
│           └── Peaks
├── scChIPseq
│   ├── MM468
│   │   ├── BigWigs
│   │   ├── Count_Matrices
│   │   └── Raw_Counts
│   │       ├── K27
│   │       └── K4
│   └── PDX
│       ├── BigWigs
│       │   ├── HBCx95_m43_UNT_H3K27me3.bw
│       │   └── HBCx95_m43_UNT_H3K4me3.bw
│       └── Count_Matrices
│           ├── HBCx95_m43_UNT_H3K27me3_TSS.tsv.gz
│           └── HBCx95_m43_UNT_H3K4me3_TSS.tsv.gz
└── scRNAseq
    ├── MM468
    └── PDX


    
```
Please refer to the scripts if you have doubts where you should place your input 
files.

## 1.0 scRNA-seq

### PDX models
Note that:  
* BC976 stands for PDX_95/Patient_95  
* BC408 stands for PDX_39/Patient_39  
* BC1224 stands for PDX_172/Patient_172  

### MDA-MB-468 model
There are 3 sub-analyses in this directory: the main analysis ('1.Persister') and
two sub analyses in response to epigenomic drug treatment ('2.UNC' and '3.UNC_5FU').
First run the 1.0 and 2.0 QC scripts.
Then you can run the script of each analyses indepentently by respecting the order
in each given analyses. Note that you must first run the single-cell ChIPseq analyse
before running the script '*5.0_ComparisonChIPseq.R*'.


## 2.0 scChIP-seq
Run the scripts of H3K27me3 and H3K4me3 in any order.

## 3.0 ChIPreChIP (a.k.a. Sequential ChIP-seq)
Run these scripts to retrieve bivalent genes as well as bivalent pathways.

## 4.0 bulk ChIPseq
These scripts are mainly to produce snapshots of specific genes from the bigwigs.


## Additional Files:

- input/bulk_ChIPseq/MM468/chromatin_indexing_qc.csv -> ratios of Chromatin Indexing IP / input
- input/scChIPseq/MM468/Raw_Counts/ -> raw counts used to calculate FrIP for scChIP
- input/scChIPseq/MM468/BigWigs/ clusters C2 / C3 / C4 to produce Extened Figure 5h  
