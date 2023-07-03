# polio_econ
This repository includes code for reproducing the results presented in the manuscript "Outbreak risks, cases, and costs of vaccination strategies against wild poliomyelitis in polio-free settings: a modelling study" by M. Auzenbergs et al. 2023

To run the code on your own, the file "psia_model.R" can be used. If you wish to achieve a steady state before introducing the first infection, run "steady_state.R" first. Once all desired models have been run, model outputs can be processed using "processing_costs.R". Alternatively, the .csv files can be read into "processing_costs.R" to validate and reproduce economic results.

Sample data from three vaccination strategies has been included in this respository: (I) RI + oSIAs, (II) RI + oSIAs + annual pSIAs, (III) RI + oSIAs + biannual pSIAs. The sample data includes model outputs with baseline RI coverage levels ranging form 25-100% and each model was run for 1,000 simulations. 

The data has undergone some post model processing to include the total number of outbreaks and expected VAPP cases. These assumptions are clarified in the manuscript.
