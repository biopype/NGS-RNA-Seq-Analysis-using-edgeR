# Differential Gene Expression Analysis using edgeR  
## Description
Crohn’s Disease (CD) is a chronic inflammatory disorder affecting the digestive tract.  
This project analyzes **RNA-Seq data from ileal biopsies** of pediatric CD patients and healthy controls.

- Reads were aligned to the human reference genome (Ensembl release 95) using HISAT2  
- Gene-level counts were generated using FeatureCounts  


## Objective
To perform **Differential Gene Expression (DGE)** analysis using **R and edgeR**, and interpret biologically significant patterns between:
- Crohn’s Disease (CD) samples  
- Control samples  


## Libraries Used

| Library | Purpose |
|--------|--------|
| edgeR | Core differential expression analysis |
| limma | Linear modeling support |
| EnhancedVolcano | Volcano plot visualization |
| pheatmap | Heatmap generation |
| ggplot2 | General data visualization |
| dplyr | Data manipulation and cleaning |
| RColorBrewer | Color palettes for plots |


## Dataset Overview

- **Total genes:** 65,217  
- **Total samples:** 304  
- **Metadata column:** `Category` (CD vs Control)  

Preprocessing: 
- Set **control as baseline/reference**  


## Data Analysis Workflow

### 1. Preprocessing
- Filtering low-expression genes  
- TMM normalization (library size correction)  
- Log transformation → **logCPM (Counts Per Million)**  


### 2. Sample-to-Sample Analysis
- Distance heatmap generated  
- Too dense due to 304 samples → excluded  


### 3. MDS Plot (Outlier Detection)

MDS (Multi-Dimensional Scaling) was used to visualize clustering.

**Observations:**
- Control samples form a **tight cluster**  
- CD samples show **higher variability**  

Potential outliers identified beyond **x > 2.5**

#### Figure 1: MDS Plot
<img width="1482" height="912" alt="MDS plot without labels" src="https://github.com/user-attachments/assets/4cbf5411-bda4-4955-b0bf-71d3ad4c90bb" />

#### Figure 2: Outlier-labeled MDS Plot
<img width="1482" height="912" alt="plot" src="https://github.com/user-attachments/assets/3e209dad-7fe3-4a4c-8fe7-a342c5dcf3af" />

**Outliers:**
- SRR5862089  
- SRR5862070  
- SRR5862173  
- SRR5862086  
- SRR5862219  
- SRR5862092  

Outliers were **not removed**


## Differential Gene Expression Analysis

### Method
- Quasi-Likelihood (QL) model  
- QL F-test  
- FDR correction  


## Results

### Volcano Plots

#### Figure 3: Volcano Plot
<img width="1902" height="967" alt="volcano without labels" src="https://github.com/user-attachments/assets/1a691baf-41f5-4204-a4cc-11a12135978e" />

### Key Findings

Thresholds:
- **p-value < 1e-20**
- **log₂FC > 5.5**

Results:
-  7 upregulated genes  
-  0 downregulated genes  


### Cutoff Comparison

| p-value | log₂FC | Upregulated | Downregulated |
|--------|--------|------------|--------------|
| 1e-15 | 5.0 | 18 | 0 |
| 1e-20 | 5.0 | 10 | 0 |
| 1e-20 | 5.5 | 7 | 0 |


## Heatmap of Significant Genes

Top 30 DEGs (by FDR):

#### Figure 5: Heatmap
<img width="1920" height="1013" alt="top 30 DEGs" src="https://github.com/user-attachments/assets/a93bfba5-3fc9-453d-af70-d7d09a96949d" />

**Insights:**
- Clear expression differences between groups  
- Strong clustering patterns  


## Conclusion

- Identified **significant upregulated genes** in CD  
- No downregulated genes under strict thresholds  
- CD samples show **higher biological variability**  
- Visualizations (MDS, Volcano, Heatmap) helped interpret results  


## Future Work
- GO / KEGG pathway analysis  
- Functional annotation  
- Machine learning models  
- Biomarker discovery  


## Notes
- Analysis performed in **R (edgeR pipeline)**  
- Fully reproducible workflow  
