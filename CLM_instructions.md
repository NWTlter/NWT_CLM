# Creating a Case for Niwot BGC simulations
## Get CLM (with proper code changes)
# Create surface datasets [I've already done this]
## 1) Creating the case
Move into the code base directory

```bash
cd /glade/u/home/hholland/clm5.0-master/cime/scripts
```
Create a new case
```bash
./create_newcase --compset  2000_DATM%GSWP3v1_CLM50%BGC_SICE_SOCN_SROF_SGLC_SWAV --res f09_g17 --case /glade/u/home/$USER/clm_point_sim_cases/2000datm_CLM50bgc_nwt_DM --run-unsupported --project P93300041
```
## 2) Set up the datm streams

There should be three datm streams but they will all point to the same files which will include all of the forcing variables.
The three files are as follows:
```bash
user_datm.streams.txt.CLMGSWP3v1.Precip
user_datm.streams.txt.CLMGSWP3v1.Solar
user_datm.streams.txt.CLMGSWP3v1.TPQW
```

In each datm stream several things need to be changed. 
1) Verify that the location of the domain file is correct
2) Verify that the location of the forcing files are correct
3) Verify that the file names match the forcing files

As an example here is what the Precip user_datm file looks like:
```xml
<?xml version="1.0"?>
<file id="stream" version="1.0">
<dataSource>
   GENERIC
</dataSource>
<domainInfo>
  <variableNames>
     time    time
        xc      lon
        yc      lat
        area    area
        mask    mask
  </variableNames>
  <filePath>
     /glade/work/hholland/single_point_nwt_neon_forcings
  </filePath>
  <fileNames>
     domain.lnd.fv0.9x1.25_gx1v7_254.42_40.05.151020.nc
  </fileNames>
</domainInfo>
<fieldInfo>
   <variableNames>
     PRECTmms precn
   </variableNames>
   <filePath>
     /glade/work/hholland/single_point_nwt_neon_forcings/datmdata_old_2008-2013/out_2008_2013_posFSDS
   </filePath>
   <fileNames>
2008-01.nc
2008-02.nc
2008-03.nc
2008-04.nc
2008-05.nc
2008-06.nc
2008-07.nc
2008-08.nc
2008-09.nc
2008-10.nc
2008-11.nc
2008-12.nc
2009-01.nc
2009-02.nc
2009-03.nc
.... (list as many forcing files as you have here)
2013-10.nc
2013-11.nc
2013-12.nc
   </fileNames>
   <offset>
      0
   </offset>
</fieldInfo>
</file>

```

Copy the parameter file into the case directory:
```bash
cp ../Surface_and_parameter_data/clm5_params.c200519_modfrootleaf_modleafcn.nc .
```
## 3) Setup the case computational requirements

Modify the *.xml files to specify the run conditions and computer-use requirements
```bash
./xmlchange --file env_mach_pes.xml --id COST_PES --val 36 # cores used relative to # tasks
./xmlchange --file env_mach_pes.xml --id TOTALPES --val 1 # total physical cores used
./xmlchange --file env_mach_pes.xml --id NTASKS --val 1 # number of mpi tasks
./xmlchange --file env_mach_pes.xml --id NTASKS_PER_INST --val 1 # ntasks per instance of each component
./xmlchange --file env_mach_pes.xml --id ROOTPE --val 0 # global mpi task of component

# Changes to all xml files
./xmlchange MPILIB=mpi-serial # in any xml file with MPILIB, change value to mpi-serial

# Changes to env_workflow.xml
./xmlchange --file env_workflow.xml --id JOB_WALLCLOCK_TIME --val 6:00:00 # time alloted for job to run - much longer for BGC run.
./xmlchange --file env_workflow.xml --id JOB_QUEUE --val share # run on shared queue - this will run faster because you can share nodes with other users

./case.setup

```
## 4) Setup the case basic simulation conditions
Modify the *.xml files to specify the simulation conditions
```bash
# Starting conditions
./xmlchange --file env_run.xml  --id CLM_FORCE_COLDSTART --val on # CLM forced to do a cold start with arbitrary initial conditions (also used for spin-up). finidat is set to blanks
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_START --val 2008 # starting year to loop data over in spinup - should correspond to first year of forcing data
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_END --val 2016 # endingyear to loop data over in spinup - should correspond to last year of forcing data
./xmlchange --file env_run.xml --id STOP_OPTION --val nyears # the run length; works with STOP_N and STOP_DATE

# Set conditions for BCG run for Accelerated Decomposition during startup
./xmlchange --file env_run.xml --id CLM_ACCELERATED_SPINUP --val on # use the accelerated decomposition mode for spinup
./xmlchange --file env_run.xml --id STOP_N --val 50 # gives the numerical value of STOP_OPTION; if STOP_OPTION is nyears and STOP_N is 5, the model runs for 5 years. 
./xmlchange --file env_run.xml --id REST_N --val 50 # sets the model restart rights; takes the value of REST_OPTION, which has the same units as the STOP_N by default. (i.e. years)
./xmlchange --file env_run.xml --id RUN_REFDATE --val 0000-01-01 # reference date for hybrid or branch runs
./xmlchange --file env_run.xml --id RUN_STARTDATE --val 0000-01-01 # the startdate - only used for startup or hybrid runs; branch runs ignore this value in favor of a reference date specified in the restart dataset. 
./xmlchange --file env_run.xml --id CONTINUE_RUN --val FALSE
./xmlchange --file env_run.xml --id RESUBMIT --val 7 # run in 8 batchs of 50 years each; 400 years total

# Finally we set the same options for the locations of the land and atmosphere domain files
# Locations of the Atmosphere and land data
# Locations of the Atmosphere and land data
./xmlchange --file env_run.xml --id ATM_DOMAIN_FILE --val domain.lnd.fv0.9x1.25_gx1v7_254.42_40.05.151020.nc
./xmlchange --file env_run.xml --id ATM_DOMAIN_PATH --val /glade/work/hholland/single_point_nwt_neon_forcings
./xmlchange --file env_run.xml --id LND_DOMAIN_FILE --val domain.lnd.fv0.9x1.25_gx1v7_254.42_40.05.151020.nc
./xmlchange --file env_run.xml --id LND_DOMAIN_PATH --val /glade/work/hholland/single_point_nwt_neon_forcings

```

## 4) Change the user namelist files
First, set the surface dataset file. This file differs for the DM, MM, WM, and SB communities. 


- Several changes have been made to the surface datasets. All have had their clay/sand percentages adjusted to mimic the rocky soil more accurately. 
- In addition, a different soil depth must be used for each vegetation community: 
(Table adapted  from Wieder et al. 2007)

| Community         | Soil Depth (cm) | Snow Conditions (% relative to obs) | Surface file name |
| ----------------- | --------------- | ----------------------------------- | ----------------- | 
| Fellfield (FF)    | 70              | 10, but 25 March, April, May (MAM)  | `surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_70cm_soildepth.nc`  |
| Dry meadow (DM)   | 100             | 10, but 25 March, April, May (MAM)  |`surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_100cm_soildepth.nc` |
| Moist meadow (MM) | 130             | 100                                 | `surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_130cm_soildepth.nc` |
| Wet meadow (WM)   | 100             | 75, + runoff from simulation of MM  |`surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_100cm_soildepth.nc` |
| Snowbed (SB)      | 70              | 200                                 | `surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_70cm_soildepth.nc`  |


Use the following code to set the surface dataset file, changing out the location of the surface dataset depending on the vegetation community you are using.
```bash
echo "fsurdat= '/glade/u/home/hholland/clm_point_sim_cases/Surface_and_parameter_data/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_100cm_soildepth.nc'" >> user_nl_clm
```

Second, set the location of the modified parameter file. This file contains important changes to leafcn and frootleaf; it must be in the case directory
```bash
echo "paramfile= '/glade/u/home/hholland/clm_point_sim_cases/2000datm_CLM50bgc_nwt_DM/clm5_params.c200519_modfrootleaf_modleafcn.nc'" >> user_nl_clm
```

Next, set the user_nl_datm files so that the atmosphere cycles and uses nearest neighbor cells for deposition.
```bash
# Tell the atmosphere that when in doubt it should use the closest neighbor cell to our location for N deposition/aerosols etc.
echo "mapalgo= 'nn','nn','nn','nn','nn'" >> user_nl_datm
echo "taxmode= 'cycle','cycle','cycle','cycle','cycle'" >> user_nl_datm
```


Finally, set the output variables and frequency in user_nl_clm
```bash
echo "hist_mfilt= 20" >> user_nl_clm # concatenate 20 output files into one file; 20 years per file
echo "hist_nhtfrq= -8760" >> user_nl_clm # record data every year

# These lines reduce the amount of times .h0. files are written.
echo "hist_empty_htapes= .true." >> user_nl_clm
# What history variables should be recorded from simulation
echo "hist_fincl1  ='TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC', 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO'" >> user_nl_clm
```

## 5) Run the case in Advanced Decomposition mode
```bash
./case.build
./case.submit
```

## 6) Check Spinup

Copy the spinup script into the case directory and modify its options

```bash
# location of script:
cp ~/clm5.0-master/tools/contrib/SpinupStability.ncl .
```
Modify script: 
```bash
;  SPT (single point) NWT BGC 2008-2011
  caseid = "2000datm_CLM50bgc_nwt_DM" ; CHANGE ME to case name
  username = "hholland"
  annual_hist = False
  region = "SPT"                 ; Global, Arctic, or SPT (single point)
  subper = 20                    ; Subsampling period in years
  
```

Run the script
```bash
module load ncl
ncl SpinupSustainability.ncl
```
If the spinup looks good, proceed. Otherwise, run for longer until it is spun up or troubleshoot the model if things look fishy.

## 7) Run in PostAD mode

Turn off advanced decomposition mode for 200 years
```bash
./xmlchange --file env_run.xml --id CLM_FORCE_COLDSTART --val off # was set to: on; start from proscribed conditions set by the end of the AD mode run
./xmlchange --file env_run.xml --id CLM_ACCELERATED_SPINUP --val off # was set to: on; turn off AD mode

```
Set the simulation to run for 200 years starting from where it left off before
```bash
./xmlchange --file env_run.xml --id STOP_N --val 50 
./xmlchange --file env_run.xml --id REST_N --val 50 
./xmlchange --file env_run.xml --id CONTINUE_RUN --val FALSE
./xmlchange --file env_run.xml --id RESUBMIT --val 4 # run in 4 batchs of 50 years each; 200 years total
./xmlchange --file env_run.xml --id RUN_REFDATE --val 0400-01-01 # was set to: 0000-01-01
./xmlchange --file env_run.xml --id RUN_STARTDATE --val 0400-01-01 # was set to: 0000-01-01

```

./xmlchange --file env_run.xml --id RESUBMIT --val 1 # run in 4 batchs of 50 years each; 200 years total
./xmlchange --file env_run.xml --id RUN_REFDATE --val 0500-01-01 # was set to: 0000-01-01
./xmlchange --file env_run.xml --id RUN_STARTDATE --val 0500-01-01 # was set to: 0000-01-01

Tell the model where to find the proscribed conditions with which to start the model.
```bash
echo "finidat= '/glade/scratch/hholland/archive/2000datm_CLM50bgc_nwt_DM/rest/0400-01-01-00000/2000datm_CLM50bgc_nwt_DM.clm2.r.0400-01-01-00000.nc'" >> user_nl_clm
```

Resubmit the case
```bash
./case.submit
```

## 8) Check the spinup again

Run the script
```bash
module load ncl
ncl SpinupSustainability.ncl
```
If spinup looks good, proceed with transient run.

## 9) Run transient run with forcings

Point to new restart file in user_nl_clm
```bash
echo "finidat= '/glade/scratch/hholland/archive/2000datm_CLM50bgc_nwt_DM/rest/0650-01-01-00000/2000datm_CLM50bgc_nwt_DM.clm2.r.0650-01-01-00000.nc'" >> user_nl_clm
```

Change model output to half-hourly data in user_nl_clm
```bash
echo "hist_nhtfrq = 0,1" >> user_nl_clm
echo "hist_mfilt  = 1,62497" >> user_nl_clm   #writes out monthly of 42 variables (below)
# Note - Removed BTRAN since it doesn't exist in CLM5 and fixed typo in FSNO
 echo "hist_fincl2  = 'SNOWICE','FSDS','FLDS','FSR','FSA','FIRE','FIRA','FSH','FCTR','FCEV','FGEV','FGR','FGR12','FSM','TSOI','COSZEN','RAIN','SNOW','H2OSOI','WA','ZWT','ELAI','FPSN','TV','RSSUN','RSSHA','FSH_G','RHAF','RH_LEAF','RH','T10','TG','SABG','SABV','EFLX_LH_TOT','SNOW_DEPTH','SNOWLIQ','SNOWDP','QRUNOFF','INT_SNOW','FSNO'" >> user_nl_clm

```

Set env_run.xml up for transient run:
```bash
./xmlchange --file env_run.xml --id STOP_N --val 6
./xmlchange --file env_run.xml --id RUN_STARTDATE --val 2008-01-01 
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_END --val 2012

```

Check that your namelist changes work
```bash
./preview_namelists
```

# To Clone case
./create_clone --clone /glade/u/home/hholland/clm_point_sim_cases/2000datm_CLM50bgc_nwt_DM --case /glade/work/wwieder/2000datm_CLM50bgc_nwt_DMnewInupt