FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL name="jorainer/spectripyworkshop" \
      url="https://github.com/jorainer/SpectriPyWorkshop" \
      maintainer="johannes.rainer@eurac.edu" \
      description="Docker container to run the SpectriPy package demo workshop. This version bases on Bioconductor release 3.23." \
      license="Artistic-2.0"

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio --exclude=./scripts/* . /home/rstudio/

## Global installation of required packages
RUN Rscript -e "BiocManager::install(c('mzR', 'MsDataHub', 'BiocStyle', 'MsExperiment', 'xcms') , ask = FALSE, dependencies = c('Depends', 'Imports'), build_vignettes = FALSE)"

## Install keyring package from github using pak
RUN Rscript -e "install.packages('pak');pak::pak('r-lib/keyring', ask = FALSE)"

## Install package dependencies
RUN Rscript -e "pak::local_install_dev_deps(ask = FALSE)"

## Install the current package with vignettes
RUN Rscript -e "pak::local_install()"

## Install SpectriPy and caching files for rstudio user
USER rstudio
RUN Rscript -e "library(SpectriPy);library(MsDataHub);MsDataHub::PestMix1_DDA.mzML()"
RUN Rscript -e "quarto::quarto_render('vignettes/combined-r-python-ms-analysis.qmd')"

## root user needed for rstudio server properly working
USER root

## Clean up
RUN find vignettes/ -name "*.html" -type f -delete && find vignettes/ -name "*_files" -type d -exec rm -r {} + && \
    rm -rf /tmp/*
