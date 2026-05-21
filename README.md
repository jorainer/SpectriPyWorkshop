# SpectriPy: Criss-Crossing R and Python for Powerful Mass Spectrometry Data Analysis Workflows

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
![docker pulls](https://img.shields.io/docker/pulls/jorainer/spectripyworkshop)
[![Docker Image Version (latest by date)](https://img.shields.io/docker/v/jorainer/spectripyworkshop?label=docker%20image)](https://hub.docker.com/repository/docker/jorainer/spectripyworkshop)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20326951.svg)](https://doi.org/10.5281/zenodo.20326951)

## :snake: Welcome to the *SpectriPy* package demo

- :dart: show how R and Python MS data analysis methods can be combined
  in a single (interactive) workflow.
- :twisted_rightwards_arrows: combine R and Python code for:
  - chromatographic peak detection of an untargeted LC-MS/MS data file
  - extract and process MS2 spectra for detected chromatographic peaks
  - download and process public reference MS2 data
  - calculate similarity between experimental and reference MS2 spectra

---

## 🛠️ Installation

### Running through docker

The workflow file along with an R runtime environment including all required
packages and the RStudio (Posit) editor are all bundled in a *docker* container.

After installation, this docker container can be run on the computer and the
code and examples from the vignettes can be evaluated within this environment
(without the need to install any additional packages or files).

-   If you don't already have, install [docker](https://www.docker.com/). Find
    installation information [here](https://docs.docker.com/desktop/).
-   Get the [docker
    image](https://hub.docker.com/r/jorainer/spectripyworkshop) of this
    tutorial e.g. from the command line with:

```
docker pull jorainer/spectripyworkshop:latest
```

-   Start the docker container, either through the Docker Desktop, or on the
    command line with

```
docker run -e PASSWORD=bioc -p 8787:8787 jorainer/spectripyworkshop:latest
```

-   Enter [`http://localhost:8787`](http://localhost:8787) in a web browser and
    log in with username `rstudio` and password `bioc`.
-   In the RStudio server version: open any of the Quarto files in the
    *vignettes* folder and evaluate the R code blocks in that document.

### Running locally

To install all required libraries locally:

```
install.packages("remotes")
install.packages("BiocManager")
BiocManager::install("jorainer/SpectriPyWorkshop", dependencies = TRUE)
```

---

## 📌 Reproducibility & Updates

We strive for **reproducibility**. This workflows is designed to
remain stable over time, allowing you to run the vignette also in future.

-   **Major updates** will be documented here.
-   **Minor updates** can be found in the [NEWS.md](NEWS.md)

---

## 🤝 Contribution

Interested in contributing? Please check out the
[**RforMassSpectrometry Contributions Guide**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).

### 📜 Code of Conduct

We follow the
[**RforMassSpectrometry Code of Conduct**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#code-of-conduct)
to maintain an inclusive and respectful community.
