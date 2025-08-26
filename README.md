
# On prior smoothing with discrete spatial data in the context of disease mapping

This repository contains the R code to fit in NIMBLE the code to replicate and reproduce the simulation study and real-data illustration of the paper entitled *"On prior smoothing with discrete spatial data in the context of disease mapping"* [(Retegui et al., 2025)](https://doi.org/10.1177/09622802251362659).


## Table of contents

1.  [Data](#Data)
2.  [R code](#Rcode)
    1.  [Within Simulation Study](#WSS)
    2.  [Across Simulation Study](#ASS)
    3.  [RealData Illustration](#illus)
3.  [Acknowledgements](#Acknowledgements)
4.  [References](#Ref)

# Data <a name="Data"/>

This folder contains the datasets and cartography files used in the simulation studies and data illustrations presented in the work.

-   `Data_Spain_47areas.Rdata`, `Data_Spain_100areas.Rdata` and `Data_Spain_300areas.Rdata`: These datasets include counts for lung cancer and the corresponding population at risk, focusing on aggregate data for females in Spain from 2019 to 2021 divided into 47, 100 and 300 areas, respectively.

-   `Carto_Spain_47areas.Rdata`, `Carto_Spain_100areas.Rdata` and `Carto_Spain_300areas.Rdata`: These cartography files include the cartography of peninsular Spain divided into 47, 100 and 300 areas, respectively.

-   `Data_England.Rdata` and `Carto_England`: Dataset including counts for lung cancer and the corresponding population at risk, focusing on males in England in 2017 and a folder containing cartography files for England.

-   `Data_SimulationStudy_Scenario1_47areas.Rdata`, `Data_SimulationStudy_Scenario1_100areas.Rdata`, `Data_SimulationStudy_Scenario1_300areas.Rdata` and `Data_SimulationStudy_Scenario1_England.Rdata`: These datasets include the simulated data for Scenario 1 when Spain is divided into 47, 100 and 300 areas, and England, respectively. Similar file naming structure, replacing `Scenario1` with `Scenario2` or `Scenario3` can be found for Scenarios 2 and 3, respectively.

Note that data for Scenarios 4 to 6 are not included, since the spatial distribution is the same as in Scenarios 1 to 3, respectively, however, the variability of the rate values is lower.

# R code <a name="Rcode"/>

This folder contains the R code to replicate and reproduce the within prior and across priors simulation studies, as well as the data illustration described in the paper. The code is organized into three subfolders, each corresponding to a specific part of the study:

1.  Within Prior Simulation Study (available [here](https://github.com/spatialstatisticsupna/Prior_Smoothing/tree/main/R/Within_SimulationStudy)). <a name="WSS"/>

    -   This folder includes code to fit the four spatial priors discussed in Section 4.1 of the paper: iCAR, LCAR, BYM, and GP.

    -   Model-fitting scripts follow the format: `Code_WithinSimStudy_iCAR.R`, `Code_WithinSimStudy_LCAR.R`, etc.

    -   Note: the LCAR model implementation is based on [beltrán-sánchez et al (2024)](https://doi.org/10.1002/sim.10166).

    -   Results are generated using `Code_Results_Within.R`.
  

2.   Across Priors Simulation Study (available [here](https://github.com/spatialstatisticsupna/Prior_Smoothing/tree/main/R/Across_SimulationStudy)). <a name="Ass"/>

        - This folder includes both model-fitting scripts and result-generation code.
  
       - Model-fitting scripts follow the format: `Code_AcrossSimStu_iid.R`, `Code_AcrossSimStu_GP.R`, etc.
  
       - Results are generated using `Code_Results_AcrossPriors.R`.


3.  Data Illustration (available [here](https://github.com/spatialstatisticsupna/Prior_Smoothing/tree/main/R/RealData_Illustration)). <a name="illus"/>

    -   This folder contains the code used to fit the models and produce the results for the real data applications.

    -   Separate scripts are provided for the two case studies: peninsular Spain and England.

# Acknowledgements <a name="Acknowledgements"/>

The work was supported by Project PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033. Garazi Retegui is supported by PhD student scholarship from the Public University of Navarra together with Banco Santander (Ayudas Predoctorales Santander UPNA 2021-2022). Jaione Etxeberria and María Dolores Ugarte would like to acknowledge support from Project UNEDPAM/PI/PR24/05A. ![plot](https://github.com/spatialstatisticsupna/Prior_Smoothing/blob/main/micin-aei.jpg) ![plot](https://github.com/spatialstatisticsupna/Prior_Smoothing/blob/main/UNED_Pamplona_2023.jpg)

# References <a name="Ref"/>
[Retegui, G., Gelfand, A., Etxeberria, J., and Ugarte, M.D. (2025). On prior smoothing with discrete spatial data in the context of disease mapping. *Statistical Methods in Medical Research*, in press, https://doi.org/10.1177/09622802251362659](https://doi.org/10.1177/09622802251362659)


