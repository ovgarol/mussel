<!-- README for Mussel Bed Morphometry and Stress Modelling 
# SPDX-FileCopyrightText: 2026 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC-BY-4.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>
-->
# Mussel Bed Morphometry and Stress Modelling

This repository contains data and R code to analyze interannual changes in mussel and oyster bed geometry (area and perimeter) and relate them to environmental stress (air temperature and salinity) in the German Wadden Sea.

## What is included

### Main script
- `script_mussel.R`
  - Reads morphometric survey data and environmental time series
  - Builds a stress index from temperature and salinity
  - Fits a coupled area-perimeter dynamical model using a genetic algorithm
  - Exports processed data tables and publication figures

### Data folder
- `data/mussel_2014.xlsx`, `data/mussel_2019.xlsx`, `data/mussel_2020.xlsx`, `data/mussel_2021.xlsx`, `data/mussel_2022.xlsx`, `data/mussel_2023.xlsx`
  - Raw benthic bed survey files used by the script (1)
- `data/Cuxhaven_DWD!Lufttemperatur.txt`
  - Air temperature observations (Cuxhaven station) (2)
- `data/LZ_AL!Salzgehalt.txt`
  - Salinity observations (Cuxhaven Alte Liebe station) (2)
- `data/LZ_4c!Salzgehalt.txt`
  - Additional salinity station file (not used in the current script) (2)
- `data/Bertolini_2020.txt` 
  - Reference curve data similar to the observed data in the trade-off figure (3)

### Generated/derived data (written by script)
- `data/mussel_beds.csv`
  - Aggregated bed geometry and metadata by date, bank, and bank type
- `data/air_temperature.csv`
  - Daily mean air temperature time series
- `data/salinity.csv`
  - Daily mean salinity time series
- `data/stress_year.csv`
  - Annual counts of days with extreme stress conditions

### Figures (written by script)
- `fig/data_plot.pdf`
- `fig/model_plot.pdf`
- `fig/trade_off.pdf`
- `fig/stress_health.pdf`

## Model Summary

The script models two state variables:
- $A$: bed area
- $P$: bed perimeter

Environmental forcing is represented with:

$$
  \text{stress}(T,S)=\frac{1}{1+e^{(10-T)}}\cdot\frac{1}{1+e^{(S-15)}}+\frac{1}{1+e^{(20-T)}}
$$

and bounded to $[0,1]$ in the implementation.

Area and perimeter are then simulated with a coupled ODE system (Euler integration), and parameters are calibrated with the GA package to reproduce observed area and perimeter trajectories for the selected bank.

## Requirements

Install these R packages:
- readxl
- quantreg
- scales
- lubridate
- deSolve
- GA
- latex2exp

## How To Run

1. Open R in this project directory.
2. Run:

```r
source("script_mussel.R")
```

or alternatively from bash:

```bash
Rscript script_mussel.R
```

Notes:
- The workflow is configured for bank code MB and bank type Austernbank in the modeling section.
- Run the script from the project directory, and ensure the `data/` and `fig/` directories are present.

## License and Credits

- License: Software Apache 2.0 and documentation CC-BY-4.0
- Copyright: 2026 Helmholtz-Zentrum Hereon GmbH
- Contributor: Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>
- Data from:
  - (1) Landesbetrieb für Küstenschutz, Nationalpark und Meeresschutz Schleswig-Holstein used under license CC-BY-4.0
  - (2) https://www.wsv.de/ used under license DL-DE->Zero-2.0. https://www.gesetze-im-internet.de/geonutzv/ und https://www.govdata.de/dl-de/zero-2-0.
  - (3) Bertolini et al. 2020 https://doi.org/10.1016/j.aquaculture.2020.735477
