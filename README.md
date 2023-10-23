GeoMx-Workflow
===
> An app to run Nanostring GeoMx DSP (Digital Spatial Profiler) Workflows

### About
This is a Docker app built to compile a set of packages to run GeoMx Workflows [vignettes](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html) which enables comprehensive QC/analysis of the GeoMx data. This app was built on [bioconductor_docker](https://hub.docker.com/r/bioconductor/bioconductor_docker) and the sessionInfo can be found at the bottom of the page.

### Features
  * Run GeoMx Workflows handy using Docker
  * Generate QC reports with diagnostic plots
  * Analyze and annotation GeoMx data to identify biological features in dataset(s)
  * Improve research reproducibility and replicability

### Run
  * Users can pull the GeoMx-Workflow image from [Docker Hub](https://hub.docker.com/r/lootpiz/geomx-workflow)
  * The following Brain GeoMx data was downloaded from [Spatial Organ Atlas](https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/)

#### Prerequisite


#### Docker
```
$ docker run -it --rm -v /USER_DIRECTORY:/home/rstudio/analysis  lootpiz/geomx-workflow
```

#### Results
Users can find the results, i.e., PDF,TXT, and RDS files, in the following folder.
```/USER_DIRECTORY/Results```

### TO-DOs
  * Users can access and modify the parameter files
  * Create HTML report template to load diagnostic plots
  * Embed analysis modules