# Penalized function-on-function partial least-squares regression


## Routines for simulation study and real-data applications


### Description

This folder contains all the routines needed for reproducing the simulation study and real data applications described in the article: "Penalized function-on-function partial least-squares regression."

Please make sure you install our `R` package `penFoFPLS` from GitHub. You can use the following command:

```
devtools::install_github("hhroig/penFoFPLS", dependencies = TRUE)
```

Also, make sure you install all the needed packages (all available in CRAN). The routines are designed to run in parallel. All the results and graphical summaries will be generated and saved into a new folder.


### Simulation study

Using `simulation_study_pFFPLS/` as working directory, source the file: `main_simulatios_call.R`. 



### Application: gait cycle data

Using `application_gait_data_pFFPLS/` as working directory, source the file: `main_app_reps_call.R`. 



### Application: bike sharing data

Using `application_bike_data_pFFPLS/` as working directory, source the file: `main_app_reps_call.R`. 


