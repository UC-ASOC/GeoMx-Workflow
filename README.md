GeoMx-Workflow
===
> An app to run Nanostring GeoMx DSP (Digital Spatial Profiler) Workflows

## About
This is a Docker app built to compile a set of packages to run GeoMx Workflows [vignettes](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html) which enables comprehensive QC/analysis of the GeoMx data. This app was built on [bioconductor_docker](https://hub.docker.com/r/bioconductor/bioconductor_docker) and the _sessionInfo_ can be found at the bottom of the page.

## Features
  * **Run GeoMx Workflows handy using Docker**
  * **Generate QC reports with diagnostic plots**
  * **Analyze and annotate GeoMx data to identify biological features in the dataset(s)**
  * **Improve research reproducibility and replicability**

## Run
  * Users can pull the GeoMx-Workflow image from [Docker Hub](https://hub.docker.com/r/lootpiz/geomx-workflow)
  * The following Brain GeoMx data was downloaded from [Spatial Organ Atlas](https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/)
  * The following NanoString Probe Kit Configuration (PKC) file was downloaded from [Nanostring](https://nanostring.com/products/geomx-digital-spatial-profiler/geomx-dsp-configuration-files/)

#### Prerequisite
  * To run the Workflow, **three** folders under the ```/USER_DIRECTORY``` folder and **one** param file are required and the folder names are case-sensitive.
    1. **DCC** : All Digital Count Conversion files
    2. **PKC** : Probe Kit Configuration file
    3. **Annot** : [Worksheet (and sample annotation) in Excel](./example/Annotation.xlsx)
    4. **parameterSettings.R** : [Parameters to run geomx-workflow](./example/parameterSettings.R)
  
```
USER_DIRECTORY
├── Annot
│   └── Annotation.xlsx
├── DCC
│   ├── DSP-1009880000092-G-A01.dcc
│   ├── DSP-1009880000092-G-A02.dcc
│   ├── DSP-1009880000092-G-A03.dcc
│   ... skipped
│   └── DSP-1012999011009-C-F04.dcc
└── PKC
    └── Hs_R_NGS_WTA_v1.0.pkc

4 directories, 259 files
```

#### Docker
```
$ docker run -it --rm -v /USER_DIRECTORY:/home/rstudio/analysis -v /parameterSettings.R:/home/rstudio/analysis/Settings/parameterSettings.R lootpiz/geomx-workflow
```

#### Results
  * Users can find the results, i.e., PDF, TXT, and RDS files, in the ```/USER_DIRECTORY/Results``` folder.

## TO-DOs
  - [ ] Create an HTML report template to load diagnostic plots
  - [ ] Map Worksheet and metadata and generate an Excel file
  - [ ] Embed analysis modules

---
## Environment
```
R version 4.3.1 (2023-06-16)
Platform: aarch64-unknown-linux-gnu (64-bit)
Running under: Ubuntu 22.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/aarch64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/aarch64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] scales_1.2.1            reshape2_1.4.4          ggsankey_0.0.99999     
 [4] ggpubr_0.6.0            dplyr_1.1.3             RColorBrewer_1.1-3     
 [7] GeoMxWorkflows_1.6.0    GeomxTools_3.4.0        NanoStringNCTools_1.8.0
[10] ggplot2_3.4.4           S4Vectors_0.38.2        Biobase_2.60.0         
[13] BiocGenerics_0.46.0     DescTools_0.99.50      

loaded via a namespace (and not attached):
  [1] rstudioapi_0.15.0       jsonlite_1.8.7          umap_0.2.10.0          
  [4] magrittr_2.0.3          ggbeeswarm_0.7.2        farver_2.1.1           
  [7] nloptr_2.0.3            rmarkdown_2.24          zlibbioc_1.46.0        
 [10] vctrs_0.6.4             minqa_1.2.6             RCurl_1.98-1.12        
 [13] askpass_1.2.0           rstatix_0.7.2           htmltools_0.5.6        
 [16] broom_1.0.5             cellranger_1.1.0        parallelly_1.36.0      
 [19] htmlwidgets_1.6.2       plyr_1.8.9              rootSolve_1.8.2.4      
 [22] uuid_1.1-1              lifecycle_1.0.3         pkgconfig_2.0.3        
 [25] Matrix_1.6-1            R6_2.5.1                fastmap_1.1.1          
 [28] GenomeInfoDbData_1.2.10 future_1.33.0           digest_0.6.33          
 [31] Exact_3.2               numDeriv_2016.8-1.1     colorspace_2.1-0       
 [34] GGally_2.1.2            reshape_0.8.9           RSpectra_0.16-1        
 [37] progressr_0.14.0        fansi_1.0.5             abind_1.4-5            
 [40] httr_1.4.7              polyclip_1.10-6         compiler_4.3.1         
 [43] proxy_0.4-27            withr_2.5.1             backports_1.4.1        
 [46] carData_3.0-5           ggforce_0.4.1           ggsignif_0.6.4         
 [49] MASS_7.3-60             openssl_2.1.0           rjson_0.2.21           
 [52] gld_2.6.6               tools_4.3.1             vipor_0.4.5            
 [55] beeswarm_0.4.0          future.apply_1.11.0     glue_1.6.2             
 [58] nlme_3.1-163            grid_4.3.1              Rtsne_0.16             
 [61] generics_0.1.3          gtable_0.3.4            class_7.3-22           
 [64] tidyr_1.3.0             data.table_1.14.8       lmom_3.0               
 [67] sp_2.1-1                car_3.1-2               utf8_1.2.3             
 [70] XVector_0.40.0          ggrepel_0.9.4           pillar_1.9.0           
 [73] stringr_1.5.0           splines_4.3.1           tweenr_2.0.2           
 [76] lattice_0.21-8          tidyselect_1.2.0        Biostrings_2.68.1      
 [79] knitr_1.43              IRanges_2.34.1          xfun_0.40              
 [82] expm_0.999-7            pheatmap_1.0.12         stringi_1.7.12         
 [85] yaml_2.3.7              boot_1.3-28.1           evaluate_0.21          
 [88] codetools_0.2-19        tibble_3.2.1            BiocManager_1.30.22    
 [91] cli_3.6.1               reticulate_1.34.0       systemfonts_1.0.4      
 [94] munsell_0.5.0           Rcpp_1.0.11             GenomeInfoDb_1.36.4    
 [97] readxl_1.4.3            globals_0.16.2          EnvStats_2.8.1         
[100] outliers_0.15           png_0.1-8               parallel_4.3.1         
[103] bitops_1.0-7            lme4_1.1-34             listenv_0.9.0          
[106] ggthemes_4.2.4          mvtnorm_1.2-3           ggiraph_0.8.7          
[109] lmerTest_3.1-3          e1071_1.7-13            SeuratObject_4.1.4     
[112] purrr_1.0.2             crayon_1.5.2            BiocStyle_2.28.1       
[115] rlang_1.1.1             cowplot_1.1.1       
```
