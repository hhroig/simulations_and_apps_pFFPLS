---
title: "Penalized function-on-function partial least-squares regression"
subtitle: "Routines for simulation study and real-data applications"
output: html_document
date: "`r Sys.Date()`"
---



## Description

This folder contains all the routines needed for reproducing the simulation study and real data applications described in the article: "Penalized function-on-function partial least-squares regression."

Please make sure you install our `R` package `penFoFPLS` from GitHub. You can use the following command:

```
devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)
```

Also, make sure you install all the needed packages (all available in CRAN). The routines are designed to run in parallel. All the results and graphical summaries will be generated and saved into a new folder.


## Simulation study

Using `simulation_study_pFFPLS/` as working directory:

1. Source the file: `main.R`. 

2. Make sure to comment/un-comment the setting you wish to run:

```
# Setting 1:
LL <- 5 # number of basis for Y(q)
KK <- 7 # number of basis for X(p)
```

or


```{r set2, eval=FALSE}
# # Setting 2: 
LL <- 40 # number of basis for Y(q)
KK <- 40 # number of basis for X(p)
```


## Application: Canada weather data

Using `application_Canada_weather_pFFPLS/` as working directory:

1. Source `canada_weather.R` to study the performance of FFPLS vs pFFPLS in terms of cross-validated error. Settings can be changed if other scenarios are to be explored.

2. Source `canada_weather_no_reps.R` to obtain the 3D plots and 2D contour plots of the estimated coefficient functions. Settings can be changed if other scenarios are to be explored.



## Application: gait cycle data

Using `application_gait_data_pFFPLS/` as working directory:

1. Source `fda_gait.R` to study the performance of FFPLS vs pFFPLS in terms of cross-validated error. Settings can be changed if other scenarios are to be explored.

2. Source `fda_gait_no_reps.R` to obtain the 3D plots and 2D contour plots of the estimated coefficient functions. Settings can be changed if other scenarios are to be explored.




