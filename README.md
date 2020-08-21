# NWT_CLM
Workflow scripts for running point simulations of CLM at Niwot Ridge using Tvan forcing data

## Niwot scripts workflow:

1. Clean L1 tvan data using tvan_supplemental_cleaning.R
2. Use flow.lter.clm.R to generate netcdf forcings for the model.
3. Follow the instructions in CLM_instructions.md to run the model at Niwot Ridge
4. Run flow.obs.R to download and format observations for comparison with the model
5. Run flow.sim.R to format model output for comparisons with observations
6. Run Obs_sim_com_plots.R to create comparison plots between sims and model
