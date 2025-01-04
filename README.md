## Overview
**BayesianRhythm** is an R package designed to perform analysis on differential rhythmic gene expression. 
It includes tools for data normalization, hypothesis testing, and visualization, supporting multi-condition rhythmic gene expression analysis.

## Features
**Rhythmic Analysis**: Perform hypothesis testing and generate statistical results for gene expression data.
**Data Normalization**: Built-in tools to normalize gene expression counts.
**Multi-condition Support**: Analyze rhythmicity across multiple experimental conditions.
**Visualization Tools**: Generate fitted rhythmic expression curves for specific genes.

## Installation
```R
install.packages("devtools")
devtools::install_github("DrHuang123/project")

## Usage
### Example Rhythmic Analysis & Plot Rhythmic Expression
```R
library(BayesianRhythmicAnalysis)

data(ExampleData)
data <- ExampleData[["CountData"]]
time <- ExampleData[["time"]]
group <- ExampleData[["group"]]

result <- BayesianRhythmicAnalysis(
  data   = data,
  time   = time,
  group  = group,
  ncond  = 2,
  period = 24,
  verbose = TRUE
)

 PlotRhythmicGene(
   fitted_params = result$FittedParams,
   gene_id       = 166,
   period        = 24
 )


