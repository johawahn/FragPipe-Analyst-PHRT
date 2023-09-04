[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![R](https://img.shields.io/badge/R-%3E4.2-brightgreen)

# FragPipe-Analyst-PHRT

Modification of FragPipe-Analyst (WORK IN PROGRESS)

A tool for analyzing quantitative proteomics datasets for [FragPipe](https://fragpipe.nesvilab.org/).


## Features

- Differential expression analysis
- Enrichment analysis (GO/Pathways)
- Imputation (optional)
- Data visualization
  1. PCA
  2. Sample correlation
  3. Heatmaps
  4. Missing value inspection
  5. Sample coverage
  6. Protein intensity plots for selected protein(s)
  7. Imputation effect evaluation

- Report and multiple levels of exportable tables for further analysis
  - Table options
    - DE results
    - Unimputed data matrix: Original protein intensities before imputation
    - Imputed data matrix: Protein intensities after performing the selected imputation method

## Public servers

There are two server instances
- Production server is hosted at [https://fragpipe-analyst.org/](https://fragpipe-analyst.org/).
- Dev server is also hosted at [http://fragpipe-analyst.nesvilab.org/](http://fragpipe-analyst.nesvilab.org/).

## Install on your own machine

### Prerequisite
- R 4.2
- PDFlatex
  
### Multiple options
Once all the dependencies are installed, follow the steps below to run the server locally.
You can build it natively:

``` sh

# Please install git in your local computer to be able to clone the repository

# Clone the repository
git clone https://github.com/johawahn/FragPipe-Analyst-PHRT

# Move to the folder
cd FragPipe-Analyst-PHRT

# Inside R console or R studio
> library("shiny")
> runApp()
```

Or run it through Docker:

``` sh
# Clone the repository
git clone https://github.com/MonashProteomics/FragPipe-Analyst-PHRT.git

# Move to the folder
cd FragPipe-Analyst

# Build FragPipe-Analyst (Any name after -t)
docker buildx build -f Dockerfile.local -t fragpipe-analyst  --output=type=docker --platform=linux/amd64 .

# Run FragPipe-Analyst
docker run -it --platform=linux/amd64 -d -p 3838:3838 fragpipe-analyst

# Open local interface
http://localhost:3838/fragpipe-analyst/
```
