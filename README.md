## scRNA-seq-analysis-of-HCC1143
Analysis codes for single cell RNA sequencing data of breast cancer HCC1143 cell line

The R scripts here were used to generate the analysis of HCC1143 scRNA-seq data, which was described in the manuscript, Zhang et al., "Transcriptional rewiring of triple-negative breast cancer cells following treatment with paclitaxel reveals therapeutic vulnerabilities".
For any questions about the data and the analysis, please contact hmzhang.bioinfo@gmail.com.

# Data
All of the data were processed by the Cell Ranger 2.0 pipeline. The processed data used for the single cell RNA-seq secondary analysis are stored in below files:
“HCC1143_24h_Ctrl_UMI_count.txt”: gene expression profiles of 24h DMSO Control cells
“HCC1143_24h_Treated_UMI_count.txt”: gene expression profiles of 24h Paclitaxel treated cells
“HCC1143_72h_Ctrl_UMI_count.txt”: gene expression profiles of 72h DMSO Control cells
“HCC1143_72h_Treated_UMI_count.txt”: gene expression profiles of 72h Paclitaxel treated cells
“regev_lab_cell_cycle_genes.txt”: gene signature for S and G2M phases

For the first four files, please download from GEO https://www.ncbi.nlm.nih.gov/geo/.

# Scripts
The following R scripts are used to cluster and classify HCC1143 scRNA-seq data of 24h and 72h DMSO control and treated cells. The analysis involves stratification of cells by cell cycle phase, dimension reduction via PCA, visualization via t-SNE, and cluster analysis via SNN clustering.
•	scRNA_analysis_of_24h_cells.Rmd: main script to analyze the cell populations of 24h DMSO control and treated samples
•	scRNA_analysis_of_72h_cells.Rmd: main script to analyze the cell populations of 72h DMSO control and treated samples

# Reference
Transcriptional rewiring of triple-negative breast cancer cells following treatment with paclitaxel reveals therapeutic vulnerabilities, in review, 2020.

