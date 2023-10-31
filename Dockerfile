# Docker inheritance
FROM bioconductor/bioconductor_docker:3.17-R-4.3.1

# Container metadata
LABEL maintainer="Heewon Seo <Heewon.Seo@UCalgary.ca>"

# Update Ubuntu
RUN apt-get update \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Install Bioconductor packages
RUN R -e 'BiocManager::install("NanoStringNCTools")'
RUN R -e 'BiocManager::install("GeomxTools")'
RUN R -e 'BiocManager::install("GeoMxWorkflows")'

# Install CRAN packages
RUN R -e 'install.packages("stringr")'
RUN R -e 'install.packages("yaml")'
RUN R -e 'install.packages("dplyr")'
RUN R -e 'install.packages("scales")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("ggpubr")'
RUN R -e 'install.packages("reshape2")'
RUN R -e 'install.packages("DescTools")'
RUN R -e 'install.packages("RColorBrewer")'
RUN R -e 'install.packages("devtools")'

# Install packages from GitHub
RUN R -e 'devtools::install_github("davidsjoberg/ggsankey")'

# Create a working directory
RUN mkdir /home/rstudio/R
RUN mkdir /home/rstudio/analysis

# Copy a set of code
ADD R /home/rstudio/R

# Start app
CMD Rscript /home/rstudio/R/geoMxWorkflow.R
