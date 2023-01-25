# LMC Analysis Examples
#### This folder contains example analyses of simulation data from Bustard et al. 2020, which modeled the Large Magellanic Cloud, a dwarf satellite galaxy of the Milky Way, with a modified version of the FLASH magnetohydrodynamics code.

## Main material
- LMC_Part1.ipynb -- a Jupyter notebook with background on the simulations and working examples of data analysis, specifically the creation of mock observables from simulation output. The mock observables we generate include Faraday rotation measure (which is a probe of magnetic field strength), gamma-ray emission (which probes cosmic ray content), H-alpha emission, and ion column densities; products from this work have been published in 3 papers and have been included in successful computing proposals. We'll additionally make a few mock absorption line spectra comparable to that from the COS-Halos spectrograph aboard the Hubble Space Telescope; similar but more extensive analysis has been the basis for multiple undergraduate projects I've advised. 

- LMC_Part2.ipynb -- follow-up to part I. This notebook plots values of "ram pressure" across the LMC disk and tabulates them at specific points where there are observational sightlines (this part, followed by an analysis of correlations between local ram pressure and observed ion abundances forms the basis for a funded Hubble Space Telescope (HST) grant (PI: Yong Zheng, Co-I: Chad Bustard)

## Supplemental material
- In the `Extras` folder, there are Fortran (.F90) files that create the initial conditions and boundary conditions for the simulations. 
- There is also a `CoolingTables` folder with a `makeCoolingFile.py` Python script I created to generate a cooling curve from the Haardt and Madau cooling tables. 
- The `Proposal` folder contains the awarded computing proposal (PI: Chad Bustard) for this project, which gives more overview of the motivation for this work, as well as scaling and timing tests for our modified version of the FLASH code.
- `Image_Video_Files` contains supplemental images and videos I've loaded into the Jupyter notebooks. 
