# NWT_CLM
This repoistory contains scripts that are necessary for running and analyzing data from CLM point simulations at Niwot Ridge, using Tvan Forcing data. 

## Installing
First fork this repository:

Then clone your fork:

```
git clone https://github.com/$USER/NWT_CLM.git
```

## Creating a python environment
Create the environment with:

```bash
cd NWT_CLM
conda env create -f environment.yml
```

Then you should hopefully have **ctsm-py** available as an environment in JupyterHub

Then intall the utilities:

```bash
conda activate NWT_CLM
pip install -e .
```


## Niwot scripts workflow:

1. [Clean L1 data](#clean-l1-data) from tvan data using `clean_tvan_data.R`*
2. [Gap fill and generate .nc forcings](#generate-atmospheric-forcings-for-clm) with `prepare_forcings_for_clm.R` 
   - also requires `nwt_QRUNOFF.ipynb` to supplement wet meadow precipitation
3. [Run the model](#run-the-model) at Niwot Ridge by following the instructions in `CLM_instructions.md` 
4. [Download and format observations](#download-and-format-observations) for comparison with the model using `prepare_obs_for_comparison.R`
5. [Format model output](#format-model-output) for comparisons with observations using `prepare_sim_for_comparison.R` 
6. [Create comparison plots](#comparing-model-and-obs) between simulation and observations with `plot_obs_sim_comparisons.R`
7. [Additional comparison plots](#comparing-model-and-obs-II) between simulation and observations with `plotFluxes_TVan_CLM.ipynb` & `plot_NWTcommunities_CLM.ipynb`

*This script will be rendered obsolete once the Tvan data is available on AmeriFlux 

![**Figure 1**: NWT_CLM workflow overview](./images/Conceptual_diagram.png)



## How to run each script

# Clean L1 data
### 1. `clean_tvan_data.R` 
This script cleans up Tvan L1 data that has been produced with the Niwot LTER `tvan_L1_preprocess.R` script from the Niwot LTER's repository. It is a temporary script meant to add a few extra cleaning steps to the L1 tvan data, until that data can be uploaded to AmeriFlux. The script reads in ReddyProc-ready data output by the `tvan_L1_preprocess.R` script, filters several problem spots, plots yearly comparisons between the filtered and unfiltered data, downloads Saddle Met data from EDI to fill in the gaps in air temperature after 2016, and writes out the data to two files called `tvan_[tower]_[start_timestamp]_to_[end_endtimestamp]_flux_P_reddyproc_cleaned.txt`

#### Inputs
 1. ReddyProc-ready files from `tvan_L1_preprocess.R`; This script expects those files to have the following variables
  - `NEE` - Net ecosystem exchange (umolm-2s-1)
  - `LE` - Latent heat flux (Wm-2)
  - `H` - Sensible heat flux (Wm-2)
  - `Ustar` - friction velocity (ms-1)
  - `Tair` - Air temperature (degC)
  - `VPD` - Vapor pressure density (kPa)
  - `rH` - relative humidity (unitless fraction)
  - `U` - Wind speed (ms-1)
  - `P` - Atmospheric Pressure (kPa)
  - `Tsoil` - Soil temperature (degC)
  - `Year` - Year of measurement (MST)
  - `DoY` - Day of year of measurement (MST)
  - `Hour` - decimal hour of measurement (MST)
  
 2. Air temperatures taken from Saddle Met data from EDI. This is automatically downloaded from EDI. Since it is only meant to fill in the gap in air temperature at 2016, only gaps from 2016+ are filled with Saddle air temperature data. 
 
#### User Options
 - `makeplots` - should plots be made? (`TRUE` or `FALSE`)
 - `DirOutBase` - The output directory for the script. It is recommended but not required that this be the same directory that holds the Reddyproc-ready files produced by tvan_L1_preprocess.R
 - `tower` - the tower data be supplementally processed, options are "East", "West", or "Both". If "Both" the both towers will be processed at once.
 - `east_data_fp` - The location of the east tvan data filepath, use "", if `tower = "West"`
 - `west_data_fp` - The location of the west tvan data filepath, use "", if `tower = "East"`

#### Outputs
The script will create a directory called `supp_filtering` in the `DirOutBase` location and save the filtered data to that directory. A directory to hold the Saddle Met data will be created within this directory and if `makeplots = TRUE` it will also create a directory called `plots` within the `supp_filtering` directory to hold the yearly plots of each variable. 

File structure: 

```bash
<DirOutBase>
└── supp_filtering
    ├── Plots
    |   └── [variable]_yearly_plots
    ├── tvan_[tower]_[start_timestamp]_to_[end_endtimestamp]_flux_P_reddyproc_supproc.txt
    └── saddle_met_data

```

#### Script Status
This script is currently necessary to further clean the tvan forcings in preparation for `prepare_forcings_for_clm.R` in addition to gap-filling one tower with the other. After a finalized version of the tvan data is up on Ameriflux, this script may no longer be necessary or may be used only for combining the two tower's data together to fill any gaps. 

The proposed modifications for this script once the Tvan data are on Ameriflux are:  

1. Copy the `download_amflx()` function and `Handle Radiation data` section from `prepare_forcings_for_clm.R` to this script

2. Modify the copied as needed to automatically download and read-in the NR-2/NR-3 data (tvan east and west towers).

3. Verify that the rest of the script as written works with the new format of tvan data.

[top](#nwt_clm)


# Generate atmospheric forcings for CLM
### 2. `prepare_forcings_for_clm.R`
The `prepare_forcings_for_clm.R` script generates atmospheric forcings for CLM from Niwot Ridge. It assembles the forcings with observational data from three sources. Daily precipitation data from the saddle that has been distributed into half-hourly data according to the method laid out in Wieder et al. 2017, half-hourly radiation data from the NR1 AmeriFlux tower, the rest of the forcings from the Tvan towers at Niwot.

| <span> |
| :--- |
| **NOTE:** For this script to work, the user must have an AmeriFlux username and account. |
| <span> |

#### Inputs
 1. Tvan data from the Niwot Ridge Tvan towers, either tower can be used, or both. If both are used, then one tower will be used to gap-fill the other. For the variables that are fed into the model, the two towers have good congruence. Right now, the user specifies the location the data generated by `supplemental_cleaning.R`, but eventually, the data will be on Ameriflux and the `download_amflx()` function can be used to download this data. 
 
 2. Saddle daily precipitation data, these are automatically downloaded from EDI, as are C1 precipitation data from USCRN. The Saddle precip data has a blowing snow correction applied to months of Oct-May, see William et al. (1988). The C1 data are used to distribute the daily precipitation from the Saddle proportionally into 30-minute timesteps. 
 
 3. AmeriFlux NR1 tower radiation data, these are automatically downloaded and used to provide short and long-wave radiation data for the forcings. 


#### User Options
 - `makeplots` - should plots be made? (`TRUE` or `FALSE`)
 - `DirOutBase` - the base directory for output, tagged with time and date version
 - `DirDnld` - Directory to download precipitation and radidation data to
 - `getNewData` - flag to determine if a newer version of precip data be automatically downloaded if one is available.
 - `amf_usr` - AmeriFlux username; NOTE: you cannot download Ameriflux data without a valid username to create an account, visit the Ameriflux website: https://ameriflux.lbl.gov/ Please also read their data-use policy, by downloading their data you are agreeing to follow it. The policy can be found here: https://ameriflux.lbl.gov/data/data-policy/
 - `tower` - the Tvan tower that will be used for the forcings. Options are "East", "West", or "Both". If "Both" the one tower will be used to gapfill the other tower basetower provides which tower is the baseline that will be filled with the other tower. Currently the East tower record is more complete and has fewer gaps and errors, so it is being used as the basetower.
 - `basetower` - The tower that will be used as the default, any gaps, will be filled with the other tower if tower is set to "Both"; Options are: "East" or "West". Recommended default is "East"
 - `east_data_fp` - The location of the east tvan data filepath, use "", if tower = "West". This is the location of the East tower output from `supplemental_cleaning.R` 
 - `west_data_fp` - The location of the east tvan data filepath, use "", if tower = "East". This is the location of the West tower output from `supplemental_cleaning.R` 
 
Options currently under development:

- `simulated_runoff_fp` - A character string specifying the location of the simulated runoff data from a moist meadow simulation to be added to Wet meadow precipitation. 

#### Outputs

**The `data` folder**

 - Five folders containing netcdf files for each of the precipitation regimes for the 4 vegetation communities, plus a folder with netcdf files of unaltered precipitation ("original") are produced. 
 - 3 text files are also produced:
   - `tvan_forcing_data_[tower]_[start_date]_[end_date].txt` - The ungapfilled data (prior to ReddyProc gapfilling); if `tower = "Both"` then this will be the data after both towers have been combined. 
   - `tvan_forcing_data_flagged_both_towers_[start_date]_[end_date].txt` - Flags indicating where data from one tower has been used to gap-fill another tower if `tower = "Both"`
   - `tvan_forcing_data_[tower]_[start_date]_[end_date].txt` - The fully gap-filled (with ReddyProc dataset including a precipitation column for each vegetation community's modified precipitation value)
   
**The `plots` folder**

 - `yearly_gap_plots_[year].png` - yearly plots showing gaps in the data prior to gap-filling of any kind
 - `all_years_gap_plots.png` - Plots showing gaps in the data prior to gap-filling for the whole period of forcing data.
 - `[start_date]_[end_date]_required_forcings_postgapfilling.png` - plots showing gaps in data after gap-filling for the whole period of forcing data
 - `[year]_yearly_gap_plots_postgapfilling.png` - yearly plots showing any gaps in the data after gap-filling is complete. 

**Directory structure:**

```bash
<DirOutBase>
└── <data_version>
    ├── data
    │   ├── tvan_forcing_data_both_towers_2007-05-11_2020-08-11.txt
    │   ├── tvan_forcing_data_flagged_both_towers_2007-05-11_2020-08-11.txt
    │   ├── tvan_forcing_data_precip_mods_both_towers_2007-05-11_2020-08-11.txt
    │   ├── dry_meadow
    │   │   ├── 2007-05.nc
    |   |   ....
    │   │   └── 2020-08.nc
    │   ├── fell_field
    │   │   ├── 2007-05.nc
    |   |   ....
    │   │   └── 2020-08.nc
    │   ├── moist_meadow
    │   │   ├── 2007-05.nc
    |   |   ....
    │   │   └── 2020-08.nc
    │   ├── original
    │   │   ├── 2007-05.nc
    |   |   ....
    │   │   └── 2020-08.nc
    │   ├── snow_bed
    │   │   ├── 2007-05.nc
    |   |   ....
    │   │   └── 2020-08.nc
    │   └── wet_meadow
    │       ├── 2007-05.nc
    |       ....
    │       └── 2020-08.nc
    └── plots
        ├── 2007-05-11_2020-08-10_required_forcing_postgapfilling.png
        ├── 2007_yearly_gap_plots_postgapfilling.png
        ....
        ├── 2020_yearly_gap_plots_postgapfilling.png
        ├── all_years_gap_plots.png
        ├── yearly_gap_plots_2007.png
        ....
        └── yearly_gap_plots_2020.png

```

#### Script Status
This script is mostly done, but could be improved by including an option to automatically create simulated run-off for the moist meadow community. The option is currently written into the user parameters as a place-holder, but is not implemented. The script expects the run-off data to be formatted as follows:

>With this option the user would provide a data file from a simulated Moist Meadow run that
contains two columns, a timestamp column (every timestamp represents the state at the *end* of the 30 minute sampling period) called "timestamp", and a column containing the runoff amounts in mm/s from a Moist Meadow  simulation. If provided, this data will be added to the Wet meadow precipitation. If not provided, wet meadow precipitation will be 75% of observed precipitation without any added runnoff. 

Steps to implement: 

 1. Proposed but untested code for this option exists in the script at lines 1727-1733 in the `Prepare 4 different precipitation regimes for the different vegetation communities` section. It needs to be tested and verified to be working.

 2. Currently there is no code to load in the simulated runoff file, that needs to be written. Ideally this data would be loaded in at the beginning of the `Prepare 4 different precipitation regimes for the different vegetation communities` section.

 3. Another possible improvement would be the creation of a function to extract the simulated runoff from a netcdf file. Currently the code as written expects a data-frame and the production of that dataframe is left to the user.

[top](#nwt_clm)

# Run the model
### 3. `CLM_instructions`
Now it's time to run CLM... see 
After running the model we'll compare with observations. First:

- The plotting scripts expect 30 minute output from the model
- Multiple years can be contatinated together
- The large dataset can be moved to your local machine for analyses

[top](#nwt_clm)


# Download and format observations
### 4. `prepare_obs_for_comparison.R`
Workflow for collating NIWOT LTER data in preparation to compare observations to simulated data. The purpose of this script is to read in observational data from Niwot, and summarize it by three time-levels: Diurnal by season, daily (day of year), and annually.

#### Inputs
 - Saddle senor network data from EDI (EDI ID: 210)
 - Saddle grid snow-depth data from EDI (EDI ID: 31)
 - Saddle productivity data from EDI (EDI ID: 16)
 - Saddle sensor network vegetation surveys from EDI (EDI ID: 191)
 - Tvan soil moisture and soil temperature data for fell-field vegetation runs
 - Tvan fluxes - from gap-filled data that was used to create netcdfs; to compare to modelled forcings 
 
 
#### User Options
 - `ver` - the data version, be default, the current system date and times
 - `DirOutBase` - the base location of the output data
 - `DirDnld` - Directory to download observation data to
 - `getNewData` - option to download the newest version of EDI data if one is available. Options: `TRUE` or `FALSE`
 - `tvan_data_fp` - Location of the tvan data that was used to create forcing files; this should be the "precip_mods" text file generated by `prepare_forcings_for_clm.R`. It must be ReddyProc-processed data because GPP is taken from this dataset. 
 - `tvan_data_soil` - Location of tvan data with soil information; This is the "Ameriflux-ready" version of the tvan data; Eventually it will need to be replaced with downloading the AmeriFlux Tvan data, once that data has been uploaded to AmeriFlux. Note: Tvan soil temperature data probes from East tower do not work, so please give west tower tvan data location


#### Outputs
This script outputs three files, the Diurnal-seasonal, daily, and yearly summaries of the observations. 

 - Diurnal-seasonal data: averaged by hour of the day across years by season. Includes mean and standard deviations of: NetRadiation (radNet), Sensible heat flux (H), Latent heat flux (LE), gross primary production (GPP)
 - Daily Data: Daily (day-of-year) means and standard deviations of soil moisture, soil temperature, and GPP; Daily snow depth at each plot.
 - Annual data: Annual means and standard deviations of production and biomass data

#### Script Status
This script could be improved by expanding it to include observations of saddle grid productivity data and creating as summary of annual GPP/NPP/ANPP observations (as was done in Wieder et al. 2017 fig 5). There is some code for this in the `Handle Saddle Grid Productivity data` section, but it is not currently used by the downstream plotting script, and the units are likely incorrect. Determining how best to separate the available data into GPP/NPP/ANPP is also another challenge. 

There also appear to be some bugs with the Soil Moisture and Temperature readings. In particular the fell field soil moisture and temperature are likely unreliable (they come from the Tvan measurements).

Other observations that could be added but are not currently in the code at all: Growing season length, soil moisture stress, delta N limitation, biomass. Determining what long term data to use for these observations, downloading it, and converting it into a useable form, will be the next step to bring this script up to scratch.

In general, there is a lot of future work to be done on this script before we can recreate all of the plots in Wieder et al. 2017 and beyond.


[top](#nwt_clm)

# Format model output
### 5. `prepare_sim_for_comparison.R`
The goal of this script is to read in netcdf output data from a CLM model simulation and convert it into tab-delimited files that can be compared to observational data. The data produced is half-hourly but is also summarized at three levels: Diurnal by season, daily (day of year), and annually. Depending on the vegetation community a different level of soil will be used for the upper layer of soil moisture and soil temperature data. This is because Tvan soil data is collected at 10cm while the saddle sensor network soil data is collected at 5cm. 

#### Inputs
 - At least one netcdf file from a transient CLM point simulation at Niwot Ridge (following instructions in `CLM_instructions.md`); The netcdf file should be half-hourly output and at a minimum contain the following history fields:
 
 | Field Name  | Description                                                | Units    |
 | ----------- | ---------------------------------------------------------- | -------- |
 | FSH         | sensible heat flux                                         | W/m^2    |
 | T10         | temperature at 2m                                          | degC     |
 | RNET        | net radiation FSA-FIRA                                     | W/m^2    |
 | TSOI_2      | Soil temperature | 2nd layer; 4cm, 2-6cm layer)            | degC     |
 | TSOI_3      | Soil temperature (3rd layer; 9cm, 6-12cm layer)            | degC     |
 | TSOI_5      | Soil temperature (5th layer; 26cm, 20-32cm layer)          | degC     |
 | H2OSOI_2    | volumetric soil moisture (2nd layer; 4cm, 2-6cm layer)     | mm/mm    |
 | H2OSOI_3    | volumetric soil moisture (3rd layer; 9cm, 6-12cm layer)    | mm/mm    |
 | H2OSOI_5    | volumetric soil moisture ( 5th layer; 26cm, 20-32cm layer) | mm/mm    |
 | SNOW_DEPTH  | snow height of snow covered area                           | cm       |
 | EFLX_LH_TOT | total latent heat flux [+ to atm]                          | W/m^2    |
 | GPP         | Gross primary production                                   | gC/m^2/s |
 | AGNPP       | Aboveground net primary productivity                       | gC/m^2/s |
 | NPP         | Net primary production                                     | gC/m^2/s |
 
#### User Options

 - `DirOutBase` - Base directory for script output
 - `case_name` - The case name of the simulation, this is used for creating subdirectories to hold the output data
 - `DirIn` - The input directory where simulation data is located
 - `[veg_com]_ncdf_fp` - The name of the netcdf file from each vegetation community ("FF", "DM", "WM", "MM", and "SB"). If you don't have a netcdf for a particular community or don't want to specify it, set the file name to `""`.
 - `usr_var` - option to export extra variables that the user can choose to select from the netcdf files. Note, if the user wants mean annual summaries of ELAI or TOTVEGC, they must be specified here since they are not included in the default variables above. 

#### Outputs

A folder in the base output directory named after the case_name with:
 1. Diurnal-seasonal data: all chosen variables averaged by hour of the day across years by season
 2. Daily Data: Daily (day-of-year) means and standard deviations of all chosen variables.
 3. Annual data: Annual means and standard deviations of all variables
 4. Unsummarized data: all chosen variables at all timestamps during the simulation
 5. Unit definitions: Units for each of the variables that are written out.
 6. Mean_annaual_summaries_vars: mean annual summaries of several variables ("GPP", "NPP", "ET", "TOTVEGC") specified in the static workflow parameters. NOTE: if ET or TOTVEGC are not specified by the user in `usr_var` they will not be summarized since they are not included in the default variable list.  
 7. Max_elai_summary: a summary by vegetation community of the max elai per year and averaged over all years.

Note: For outputs 1-4, data from all available vegetation communities are concatenated and saved to a single file.

#### Script Status
This script is listed as "Done". 

[top](#nwt_clm)


# Comparing model and obs
### 6. `plot_obs_sim_comparisons.R` 
This is a script that creates plots comparing observation and simulation data after the style of Wieder et al. 2017

#### Input

 - The output from `prepare_sim_for_comparison.R`
 - The output from `prepare_obs_for_comparison.R`

#### User Options

 - `DirOutBase` - the output directory for the script
 - `sim_name` - the simulation name (for organizing output and naming outputs)
 - `DirSimIn` - the netcdf-named subdirectory containing the outputs from `prepare_sim_for_comparison.R`
 - `DirObsIn` - the version-named subdirectory containing the outputs from `prepare_obs_for_comparison.R`
 - `vegetation_com` - the vegetation community that is being compared. Options: "FF", "DM", "WM", "MM", "SB", NA;
 

#### Outputs

3 plots comparing the fluxes, soil moisture data, and snow-depth data from the simulation to observations. 

#### Script Status
This script is listed as "can be improved". It can primarily be improved by re-creating figures 5-7 of Wieder et al. 2017 (however most of these plots necessitate changes to `prepare_obs_for_comparison.R`):

 - A comparison of primary productivity outputs (GPP/NPP/ANPP) (Fig 5 of Wieder et al. 2017)
 - Plots of biomass, growing season length, N-limitation, and soil moisture stress (Fig 6 & 7 of Wieder et al. 2017). 

In addition to the plots mentioned above, only some of the plots created will compare across all vegetation community types. Ideally, the script would be modified so that each plot could be like the snow depth plot, where all communities are compared to each other. 

Finally, there are a number of plotting aesthetics that could be improved.  

[top](#nwt_clm)


# Comparing model and obs II
### 7a. `plotFluxes_TVan_CLM.ipynb` 
Also makes diel, seasonal, and annual plots of simulated and observed fluxes from Tvan
### 7b. `plot_NWTcommunities_CLM.ipynb` 
Plots simulations from 5 saddle communities and compares to observations from Saddle grid and Saddle Sensor Network.

Both scripts have the advantage of reading model output directly on /glade/scratch and are similar to scripts being developed for tower simulations at NEON sites 

*Dependencies:* `utils.py` has some additional utilities that are sparingly used

[top](#nwt_clm)

