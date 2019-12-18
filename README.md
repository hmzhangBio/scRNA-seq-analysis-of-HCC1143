# scRNA-seq-analysis-of-HCC1143
Analysis codes for single cell RNA sequencing data of breast cancer HCC1143 cell line

The R scripts here were used to generate the analysis of HCC1143 scRNA-seq data, which was described in the manuscript, Zhang et al., "Transcriptional rewiring of triple-negative breast cancer cells following treatment with paclitaxel reveals therapeutic vulnerabilities".

For any questions about the data and the analysis, please contact hmzhang.bioinfo@gmail.com.

## Data
All of the data were processed by the Cell Ranger 2.0 pipeline. The processed data used for the single cell RNA-seq secondary analysis are stored in below files:

“HCC1143_24h_Ctrl_UMI_count.txt”: gene expression profiles of 24h DMSO Control cells
“HCC1143_24h_Treated_UMI_count.txt”: gene expression profiles of 24h Paclitaxel treated cells
“HCC1143_72h_Ctrl_UMI_count.txt”: gene expression profiles of 72h DMSO Control cells
“HCC1143_72h_Treated_UMI_count.txt”: gene expression profiles of 72h Paclitaxel treated cells
“regev_lab_cell_cycle_genes.txt”: gene signature for S and G2M phases

For the first four files, please download from GEO https://www.ncbi.nlm.nih.gov/geo/.

## Scripts
The following R scripts are used to cluster and classify HCC1143 scRNA-seq data of 24h and 72h DMSO control and treated cells. The analysis involves stratification of cells by cell cycle phase, dimension reduction via PCA, visualization via t-SNE, and cluster analysis via SNN clustering.

•	scRNA_analysis_of_24h_cells.Rmd: main script to analyze the cell populations of 24h DMSO control and treated samples

•	scRNA_analysis_of_72h_cells.Rmd: main script to analyze the cell populations of 72h DMSO control and treated samples

## Reference
Transcriptional rewiring of triple-negative breast cancer cells following treatment with paclitaxel reveals therapeutic vulnerabilities, in review, 2020.

## SessionInfo()

R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X El Capitan 10.11.6

### Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib

###  locale:
en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

### attached base packages:
stats     graphics  grDevices utils     datasets  methods   base     

### other attached packages:
dplyr_0.8.0.1  Seurat_2.3.4  Matrix_1.2-14 cowplot_0.9.4 ggplot2_3.1.1

### loaded via a namespace (and not attached):
  [1] tsne_0.1-3          segmented_0.5-3.0   nlme_3.1-137        bitops_1.0-6        bit64_0.9-7         httr_1.4.0         
  [7] RColorBrewer_1.1-2  prabclus_2.2-6      tools_3.5.1         backports_1.1.2     irlba_2.3.3         R6_2.4.0           
 [13] rpart_4.1-13        KernSmooth_2.23-15  Hmisc_4.1-1         lazyeval_0.2.2      colorspace_1.4-1    trimcluster_0.1-2.1
 [19] nnet_7.3-12         npsurv_0.4-0        withr_2.1.2         tidyselect_0.2.5    gridExtra_2.3       bit_1.1-14         
 [25] compiler_3.5.1      htmlTable_1.12      hdf5r_1.0.0         diptest_0.75-7      caTools_1.17.1.2    scales_1.0.0       
 [31] checkmate_1.8.5     lmtest_0.9-36       DEoptimR_1.0-8      mvtnorm_1.0-8       robustbase_0.93-2   ggridges_0.5.1     
 [37] pbapply_1.4-0       dtw_1.20-1          proxy_0.4-22        stringr_1.4.0       digest_0.6.18       mixtools_1.1.0     
 [43] foreign_0.8-70      R.utils_2.8.0       base64enc_0.1-3     pkgconfig_2.0.2     htmltools_0.3.6     bibtex_0.4.2       
 [49] htmlwidgets_1.3     rlang_0.4.0         rstudioapi_0.7      jsonlite_1.6        zoo_1.8-5           ica_1.0-2          
 [55] mclust_5.4.1        gtools_3.8.1        acepack_1.4.1       R.oo_1.22.0         magrittr_1.5        modeltools_0.2-22  
 [61] Formula_1.2-3       lars_1.2            Rcpp_1.0.1          munsell_0.5.0       reticulate_1.12     ape_5.3            
 [67] R.methodsS3_1.7.1   stringi_1.4.3       gbRd_0.4-11         MASS_7.3-50         flexmix_2.3-14      gplots_3.0.1.1     
 [73] Rtsne_0.15          plyr_1.8.4          grid_3.5.1          parallel_3.5.1      gdata_2.18.0        crayon_1.3.4       
 [79] doSNOW_1.0.16       lattice_0.20-35     splines_3.5.1       SDMTools_1.1-221.1  knitr_1.20          pillar_1.3.1       
 [85] igraph_1.2.4.1      fpc_2.1-11.1        reshape2_1.4.3      codetools_0.2-15    stats4_3.5.1        glue_1.3.1         
 [91] lsei_1.2-0          metap_1.1           latticeExtra_0.6-28 data.table_1.12.2   png_0.1-7           Rdpack_0.11-0      
 [97] foreach_1.4.4       tidyr_0.8.3         gtable_0.3.0        RANN_2.6.1          purrr_0.3.2         kernlab_0.9-27     
[103] assertthat_0.2.1    class_7.3-14        survival_2.42-3     tibble_2.1.1        snow_0.4-3          iterators_1.0.10   
[109] cluster_2.0.7-1     fitdistrplus_1.0-14 ROCR_1.0-7  
