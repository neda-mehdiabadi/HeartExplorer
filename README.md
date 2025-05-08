The **code** here is for the interactive [HeartExplorer](https://www.HeartExplorer.org) portal that I developed. More details on the dataset and paper the web-portal is referring to are as follows:



# Defining the Fetal Gene Program at Single-Cell Resolution in Dilated Cardiomyopathy 
[![Publication: Circulation (IF 37.8)](https://img.shields.io/badge/Published%20in-Circulation%20(IF%2037.8)-red)](https://www.ahajournals.org/doi/full/10.1161/CIRCULATIONAHA.121.057763?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org) 

[![Interactive Portal](https://img.shields.io/badge/Explore%20Data-HeartExplorer.org-brightgreen)](https://www.HeartExplorer.org)  

## 📖 Overview
This project presents an in-depth analysis of **single-nucleus RNA-sequencing (snRNA-seq)** data to define the **fetal gene program** in **dilated cardiomyopathy (DCM)** at **single-cell resolution**. It supports the findings of the manuscript:

> **Defining the Fetal Gene Program at Single-Cell Resolution in Dilated Cardiomyopathy**  
> Published in *Circulation* (Impact Factor: 37.8; 2022)  
> [Link to publication](https://www.ahajournals.org/doi/full/10.1161/CIRCULATIONAHA.121.057763?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org)

## 📊 Dataset Summary
| Group         | Age Range       | Samples (n) | Nuclei Count |
|---------------|------------------|-------------|---------------|
| Fetal         | 19–20 weeks      | 3           | 27,136        |
| Non-Diseased  | 4–14 years       | 3           | 16,445        |
| DCM           | 5–10 years     | 4           | 31,593        |

- All samples are derived from the **left ventricle**.
- This is a large single nuclei RNA-seq dataset from the Heart Regeneration group at Murdoch Children's Research Institute. 
 
## 🔍 Key Analyses
-  **Cell Type Identification**  
  Annotate and cluster cell types within each group using unsupervised methods.
-  **Cell Composition Analysis**  
  Compare shifts in cell-type abundance across healthy, fetal, and DCM states.
-  **Pathway and Transcriptional Analysis**  
  Discover regulatory changes and transcription factor activity in DCM.
-  **Fetal Gene Program Recovery in DCM**  
  Identify gene expression signatures indicative of fetal-like reprogramming in disease.

## 🌐 Explore the Data
I have built an **interactive web portal** that permits interrogation of our dataset and, hence, increases the dataset’s accessibility and utility. To check our website, please visit [www.HeartExplorer.org](http://www.HeartExplorer.org/).
- Search by cell type or gene name
- View gene expression across conditions
- Perform pathway enrichment analysis across conditions
- Download selected data tables

## 📚 Citation
If you use this dataset or code, please cite:

```bibtex
  @article{MehdiabadiCirculation2022,
  author  = {Neda Rahmani Mehdiabadi and et al.},
  title   = {Defining the Fetal Gene Program at Single-Cell Resolution in Dilated Cardiomyopathy},
  journal = {Circulation},
  year    = {2022},
  doi     = {10.1161/CIRCULATIONAHA.121.057763},
  url     = {https://www.ahajournals.org/doi/full/10.1161/CIRCULATIONAHA.121.057763}
  }
```
