# Custom scripts from Fu et al

### Modeling photosynthetic responses to carbon dioxide concentration (*A/Ci* curve)

A Python script `modeling_photosynthesis_to_CO2_ACi.py`

- Cleaning the LI-6800 data from csv files in a folder
- Fitting the *Farquhar, von Caemmerer and Berry (FvCB) Model* with two terms `αG`  and `αS` , where `αG` is the proportion of glycolate carbon exported from the photorespiratory pathway as glycine, and `αS` the proportion exported as serine.
- Obtaining and summarizing the statistics of the fitting parameters

- Plotting the measured data with the three biochemical limitation curves



### Modeling photosynthetic responses during oxygen transients

An R script `modeling_photosynthesis_O2_transient.R`

- Fitting one-phase exponential functions to the time course data of net CO2 assimilation rate (A)
- Integrating the area between the fitted curve and the baseline of the new steady-state A

- Plotting the measured data with the fitted curves

