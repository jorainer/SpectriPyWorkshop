FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL name="jorainer/spectripyworkshop" \
      url="https://github.com/jorainer/SpectriPyWorkshop" \
      maintainer="johannes.rainer@eurac.edu" \
      description="Docker container to run the SpectriPy package demo workshop. This version bases on Bioconductor release 3.23." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio --exclude=./scripts/* . /home/rstudio/

## Global installation of required packages
RUN Rscript -e "BiocManager::install(c('mzR', 'MsDataHub', 'BiocStyle', 'Spectra') , ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)"

## Install keyring package from github using pak
RUN Rscript -e "install.packages('pak');pak::pak('r-lib/keyring', ask = FALSE)"

## Install SpectriPy and caching files for rstudio user
USER rstudio

RUN Rscript -e "BiocManager::install('SpectriPy')"
RUN Rscript -e "library(MsDataHub);MsDataHub::PestMix1_DDA.mzML()"

## Install the current package with vignettes
RUN Rscript -e "devtools::install('.')"

## root user needed for rstudio server properly working
USER root

## Clean up
RUN find vignettes/ -name "*.html" -type f -delete && find vignettes/ -name "*_files" -type d -exec rm -r {} + && \
    rm -rf /tmp/*
