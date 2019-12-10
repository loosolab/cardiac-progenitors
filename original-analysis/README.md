# Single cell RNA-seq and ATAC seq analysis of cardiac progenitor cell transition states and lineage settlement

### R session information

For requirements, please see the file `RNA-seq_sessionInfo.txt` and `ATAC-seq_sessionInfo.txt` for specific R packages and their versions.

### Data download

Data for the original analysis is approximately 500 MB in size and can be downloaded from a bucket of an S3 object store using [cloudyr/aws.s3](https://github.com/cloudyr/aws.s3):

From within R, install [cloudyr/aws.s3](https://github.com/cloudyr/aws.s3) as follows and set up the connection to our data server:

```
remotes::install_github("cloudyr/aws.s3")
Sys.setenv("AWS_S3_ENDPOINT" = "mpi-bn.mpg.de",
           "AWS_DEFAULT_REGION" = "s3")

if(!dir.exists("original_analysis/data")) {
    dir.create("original_analysis/data", recursive = TRUE)
}

if(!dir.exists("original_analysis/supplementary_data")) {
    dir.create("original_analysis/supplementary_data", recursive = TRUE)
}
```

Then, use the following code snipped to download all files:

```
data_files <- c("Isl1_cardiac_branch_expression.txt", "Isl1_endo_branch_expression.txt", "dev_cluster_specific_2.Rda", "dev_clustered_cluster_specific_2.Rda", "scData_filtered_2.Rda", "variability_cluster_specific_2.Rda", "variability_clustered_cluster_specific_2.Rda")

for(file in data_files) {
    aws.s3::save_object(paste0("original-data/", file),
                        "data-cpc-2018",
                        paste0("original_analysis/data/", file))
}

supplementary_data_files <- c("Isl1-DM.Rdata", "Isl1-diffExprs.Rdata", "Isl1-markers.Rdata", "Nkx2-5-DM.Rdata", "Nkx2-5-diffExprs.Rdata", "Nkx2-5-markers.Rdata", "c1_subset.Rdata", "wu.Rdata")

for(file in supplementary_data_files) {
    aws.s3::save_object(paste0("original-data/", file),
                        "data-cpc-2018",
                        paste0("original_analysis/supplementary_data/", file))
}
```

or use `aws.s3::save_object()` with the following parameters to download individual files:

Object | Target | Call
------ | ------ | ----
`Isl1_cardiac_branch_expression.txt` | `data` | `aws.s3::save_object("original-data/Isl1_cardiac_branch_expression.txt", "data-cpc-2018", "original_analysis/data/Isl1_cardiac_branch_expression.txt")`
`Isl1_endo_branch_expression.txt` | `data` | `aws.s3::save_object("original-data/Isl1_endo_branch_expression.txt", "data-cpc-2018", "original_analysis/data/Isl1_endo_branch_expression.txt")`
`dev_cluster_specific_2.Rda` | `data` | `aws.s3::save_object("original-data/dev_cluster_specific_2.Rda", "data-cpc-2018", "original_analysis/data/dev_cluster_specific_2.Rda")`
`dev_clustered_cluster_specific_2.Rda` | `data` | `aws.s3::save_object("original-data/dev_clustered_cluster_specific_2.Rda", "data-cpc-2018", "original_analysis/data/dev_clustered_cluster_specific_2.Rda")`
`scData_filtered_2.Rda` | `data` | `aws.s3::save_object("original-data/scData_filtered_2.Rda", "data-cpc-2018", "original_analysis/data/scData_filtered_2.Rda")`
`variability_cluster_specific_2.Rda` | `data` | `aws.s3::save_object("original-data/variability_cluster_specific_2.Rda", "data-cpc-2018", "original_analysis/data/variability_cluster_specific_2.Rda")`
`variability_clustered_cluster_specific_2.Rda` | `data` | `aws.s3::save_object("original-data/variability_clustered_cluster_specific_2.Rda", "data-cpc-2018", "original_analysis/data/variability_clustered_cluster_specific_2.Rda")`
`Isl1-DM.Rdata` | `supplementary_data` | `aws.s3::save_object("original-data/Isl1-DM.Rdata", "data-cpc-2018", "original_analysis/supplementary_data/Isl1-DM.Rdata")`
`Isl1-diffExprs.Rdata` | `supplementary_data` | `aws.s3::save_object("original-data/Isl1-diffExprs.Rdata", "data-cpc-2018", "original_analysis/supplementary_data/variability_clustered_cluster_specific_2.Rda")`
`Isl1-markers.Rdata` | `supplementary_data` | `aws.s3::save_object("original-data/Isl1-markers.Rdata", "data-cpc-2018", "original_analysis/supplementary_data/Isl1-markers.Rdata")`
`Isl1_cardiac_branch_expression.txt` | `supplementary_data` | `aws.s3::save_object("original-data/Isl1_cardiac_branch_expression.txt", "data-cpc-2018", "original_analysis/supplementary_data/Isl1_cardiac_branch_expression.txt")`
`Isl1_endo_branch_expression.txt` | `supplementary_data` | `aws.s3::save_object("original-data/Isl1_endo_branch_expression.txt", "data-cpc-2018", "original_analysis/supplementary_data/Isl1_endo_branch_expression.txt")`
`Nkx2-5-DM.Rdata` | `supplementary_data` | `aws.s3::save_object("original-data/Nkx2-5-DM.Rdata", "data-cpc-2018", "original_analysis/supplementary_data/Nkx2-5-DM.Rdata")`
`Nkx2-5-diffExprs.Rdata` | `supplementary_data` | `aws.s3::save_object("original-data/Nkx2-5-diffExprs.Rdata", "data-cpc-2018", "original_analysis/supplementary_data/Nkx2-5-diffExprs.Rdata")`
`Nkx2-5-markers.Rdata` | `supplementary_data` | `aws.s3::save_object("original-data/Nkx2-5-markers.Rdata", "data-cpc-2018", "original_analysis/supplementary_data/Nkx2-5-markers.Rdata")`
`c1_subset.Rdata` | `supplementary_data` | `aws.s3::save_object("original-data/c1_subset.Rdata", "data-cpc-2018", "original_analysis/supplementary_data/c1_subset.Rdata")`
`wu.Rdata` | `supplementary_data` | `aws.s3::save_object("original-data/wu.Rdata", "data-cpc-2018", "original_analysis/supplementary_data/wu.Rdata")`

Finally, set the working directory to the path `original-analysis` and continue below:

```
setwd("original-analysis")
```

### RNA-seq analysis

#### Preprocessing

R code to preprocess single-cell RNA-seq data from scratch can be found in the `src` directory. Steps (e.g. *quality control*, *filtering*, *normalization*, etc.) should be run in the order indicated by the leading counter in the filename (e.g. `00_filtering.R` before `01_normalisation.R`)

#### Figures

We additionally provide scripts to reproduce most of the main figures. The R code can be found in the `figure_src` directory and is split into individual files. Figures can be recreated without the need of *preprocessing* the data.

### ATAC-seq analysis

#### Figures

We provide scripts to reproduce figures from the ATAC-seq data analysis of the manuscript. The R code can be found in the `figure_src` directory.

### Citation

Please refer to the following research article when using data from this repository:

Guangshuai Jia, Jens Preussner, Xi Chen, Stefan Guenther, Xuejun Yuan, Michail Yekelchyk, Carsten Kuenne, Mario Looso, Yonggang Zhou, Sarah Teichmann and Thomas Braun. Single cell RNA-seq and ATAC seq analysis of cardiac progenitor cell transition states and lineage settlement. **Nature Communications** 9, 4877 (*2018*), doi: [10.1038/s41467-018-07307-6](https://doi.org/10.1038/s41467-018-07307-6).
