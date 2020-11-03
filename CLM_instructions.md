# Creating a Case for Niwot BGC simulations
## Getting and modifying the CLM code base

### CLM code changes
To run the model at Niwot Ridge, you will need to make sure that your version of CLM has several code changes implemented. The full details of these file changes can be found [here](https://github.com/hhollandmoritz/CTSM/pull/1/files)

1) **Growing season adjustment:** Phenology must be adjusted so that: 
- The growing season can start after the summer solstice, and 
- Not before soils are snow free.  Here we'll just make the min_temp for growing degree days 4C. 

(Modifications to `CNPhenologyMod.F90`)

At line 855:  
- But also see Leah's arctic phenology changes for CTSM5.1, which as a snow free constraint on decid trees.
```fortran
               ! Test to turn off growing degree-day sum, if on.
               ! This test resets the growing degree day sum if it gets past
               ! the summer solstice without reaching the threshold value.
               ! In that case, it will take until the next winter solstice
               ! before the growing degree-day summation starts again.

!WW turned off to onset can occur after summer solstice
!               if (onset_gddflag(p) == 1._r8 .and. ws_flag == 0._r8) then
!                  onset_gddflag(p) = 0._r8
!                  onset_gdd(p) = 0._r8
!               end if


               ! if the gdd flag is set, and if the soil is above freezing
               ! then accumulate growing degree days for onset trigger

               soilt = t_soisno(c, phenology_soil_layer)
               ! WW modified to increase the critical temperature only occurs after soils are unfrozen (here 4C)
               if (onset_gddflag(p) == 1.0_r8 .and. soilt > 277._r8) then !SHR_CONST_TKFRZ) then
                  onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
               end if

```

At line 921, added a note to modify `crit_dayl` for offset:

- Also see Leah's arctic phenology changes

```fortran
           ! only begin to test for offset daylength once past the summer sol
	       ! WW to correct offset, critical day length = 10.9 hours on parameter file (check for CLM5)
```

2) **Relative humidity and Photosynthesis:** Photosynthesis must be modified so that the relative humidity no longer controls photosynthesis. (Modifications to `PhotosynthesisMod.F90`)

At line 1642:

```fortran
               !now the constraint is no longer needed, Jinyun Tang
               ceair = min( eair(p),  esat_tv(p) )
               ! WW to remove rh effect on psn 
               ceair = esat_tv(p)

               if (      stomatalcond_mtd == stomatalcond_mtd_bb1987 )then
                  rh_can = ceair / esat_tv(p)
               else if ( stomatalcond_mtd == stomatalcond_mtd_medlyn2011 )then
                  ! Put some constraints on RH in the canopy when Medlyn stomatal conductance is being used
                  rh_can = max((esat_tv(p) - ceair), 50._r8) * 0.001_r8
               end if

               ! Electron transport rate for C3 plants. Convert par from W/m2 to
               ! umol photons/m**2/s using the factor 4.6
```


3) **Write out root fraction by soil layer:** modify default history fields. (Modifications to `SoilStateType.F90`)

At line 233: 

```fortran
! WW turned off to make ROOTFR = active
!    if (use_dynroot) then
       this%rootfr_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='ROOTFR', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='fraction of roots in each soil layer', &
            ptr_patch=this%rootfr_patch, default='active')
!    end if
```


At line 250:
```fortran

    if (use_cn) then
       this%rootr_patch(begp:endp,:) = spval
       call hist_addfld2d (fname='ROOTR', units='proportion', type2d='levgrnd', &
            avgflag='A', long_name='effective fraction of roots in each soil layer (SMS method)', &
            ptr_patch=this%rootr_patch, l2g_scale_type='veg', default='active') ! WW made active
    end if
```

## Generating Single Point Domain and Surface data

If you don't already have the domain and surface datasets, you will need to create them by modifying the singlept script in `clm5.0-master/tools/contrib`

1) Set the longitude and latitude to match the Tvan tower locations at Niwot Ridge.

| <span> |
| :--- |
| **NOTE:** The latitude and longitude values in this script and your netcdf forcing files *must* match! |
| <span> |

```python
plon = 254.42 # NWT LTER TVAN tower location
plat = 40.05 #  

```

2) Set the dominant PFT (plant functional type) to 12 (arctic grass)

```python
dominant_pft         = 12 # C3 arctic grass (11 is broadleaf decid. shrub - boreal)
```

3) Optional: Create alternate datm and set output directory

You can choose to create supplemental datm files for comparison from the default CLM forcing set.

```python
create_datm     = True 
```

You can choose to change the output directory to a custom location.

```python
dir_output='/glade/work/'+myname+'/single_point_nwt_forcings/'
```

4) Run the single_pt script to generate the files:

```bash
module load python
python singlept
```


## Changes to surface datasets files:

If you are not using surface datasets that are already customized for Niwot, you will need to modify them yourself:

The `singlept` script above will generate a generic surface dataset file for Niwot. This is the file that you will need to modify. 

1) Change the percent sand and clay to 39% and 23%, respectively. Use variables `PCT_SAND` and `PCT_CLAY` to do this.

```bash
module load nco # load the nco command module to access ncap2 program
# Change the percent sand and clay variables
ncap2 -s 'PCT_SAND=PCT_SAND*0+39;PCT_CLAY=PCT_CLAY*0+23' /path/to/your/surfacedata/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706.nc /path/to/your/surfacedata/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay.nc
```

2) Depending on the vegetation community that you are trying to mimic you will need to modify the depth of the soil (`zbedrock`). `zbedrock` is in units of meters. 

```bash
# Dry meadow & Wet meadow
# modify a surface dataset with proper clay/sand ratio to have bedrock at 100cm (1 m) deep
ncap2 -s 'zbedrock=zbedrock*0+1' /path/to/your/surfacedata/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay.nc /path/to/your/surfacedata/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_100cm_soildepth.nc

# Moist meadow
# 130 cm
ncap2 -s 'zbedrock=zbedrock*0+1.3' /path/to/your/surfacedata/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay.nc /path/to/your/surfacedata/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_130cm_soildepth.nc

# Snow field and Fellfield
# 70 cm
ncap2 -s 'zbedrock=zbedrock*0+0.7' /path/to/your/surfacedata/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay.nc /path/to/your/surfacedata/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_70cm_soildepth.nc

```

## 1) Creating the case
Move into the code base directory

```bash
cd /glade/u/home/hholland/clm5.0-master/cime/scripts
```
Create a new case, replace `/path/to/your/case/case_dir/` with the path to and name of the case you are running. For example:
`/glade/u/home/$USER/clm5_cases/informative_case_name_here`

```bash
./create_newcase --compset  2000_DATM%GSWP3v1_CLM50%BGC_SICE_SOCN_SROF_SGLC_SWAV --res f09_g17 --case /path/to/your/case/case_dir/ --run-unsupported --project <project_id_number>
```

Move into your case directory
```bash
cd /glade/u/home/$USER/clm5_cases/informative_case_name_here
```
## 2) Setup the case computational requirements
Modify the *.xml files to specify the run conditions and computer-use requirements
```bash
# PES changes to run light and fast!
# The we can make the following xml changes BEFORE you set up your case

./xmlchange ROOTPE=0          
./xmlchange JOB_QUEUE=share   
./xmlchange JOB_WALLCLOCK_TIME=6:00:00 
./xmlchange MPILIB=mpi-serial 
./case.setup

```

## 3) Setup the case basic simulation conditions
Modify the `env_run.xml` file to specify the simulation conditions with the following code blocks

1. Point to domain files
2. Cycle over 5 years of input data
3. Start from coldstart, in AD mode for 400 years total
- setting `RUN_REFDATE` to 0008 should let us start with 2008 forcing data?
4. then make changes to user_nl_* files, like pointing to the right surface dataset and limiting the history file output.

```bash
./xmlchange ATM_DOMAIN_FILE=domain.lnd.fv0.9x1.25_gx1v7_254.42_40.05.151020.nc
./xmlchange ATM_DOMAIN_PATH=/glade/p/cgd/tss/people/wwieder/inputdata/single_point
./xmlchange LND_DOMAIN_FILE=domain.lnd.fv0.9x1.25_gx1v7_254.42_40.05.151020.nc
./xmlchange LND_DOMAIN_PATH=/glade/p/cgd/tss/people/wwieder/inputdata/single_point

./xmlchange DATM_CLMNCEP_YR_START=2008
./xmlchange DATM_CLMNCEP_YR_END=2012
./xmlchange STOP_OPTION=nyears

./xmlchange CLM_FORCE_COLDSTART=on
./xmlchange CLM_ACCELERATED_SPINUP=on
./xmlchange STOP_N=200
./xmlchange REST_N=100
./xmlchange RUN_REFDATE=0008-01-01
./xmlchange RUN_STARTDATE=0008-01-01
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=1

echo "fsurdat = '/glade/p/cgd/tss/people/wwieder/inputdata/single_point/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_NWT_70cm_soildepth.nc'" >> user_nl_clm
echo "mapalgo = 'nn','nn','nn','nn','nn'" >> user_nl_datm
echo "hist_mfilt = 20" >> user_nl_clm
echo "hist_nhtfrq = -8760" >> user_nl_clm
echo "hist_empty_htapes = .true." >> user_nl_clm
echo "hist_fincl1 = 'TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC', 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO'" >> user_nl_clm

./preview_namelists
qcmd -- ./case.build
```

## 4) Set up the datm streams (forcing files)

If you have not already done so, move the netcdf forcing files generated by `flow.lter.clm.R` to the server. There are a different set of forcing files for each vegetation community.

In your case dirctory, there will be three datm streams but they will all point to the same netcdf files each of which includes all of the forcing variables.

The three files are as follows:
```bash
user_datm.streams.txt.CLMGSWP3v1.Precip
user_datm.streams.txt.CLMGSWP3v1.Solar
user_datm.streams.txt.CLMGSWP3v1.TPQW
```

In each datm stream several things need to be changed. 

 - [ ] Verify that the location of the domain file is correct

 - [ ] Verify that the location of the forcing files are correct

 - [ ] Verify that the file names match the forcing files

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
     /path/to/your/domain/file
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
     /path/to/your/tvan-derived/forcing/files
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

The only differences between the three user_datm.streams* files are in the variable names section. Here are those sections for the Solar and TPQW files.

Solar forcings file:

```xml
<?xml version="1.0"?>
<file id="stream" version="1.0">
....
<fieldInfo>
   <variableNames>
     FSDS swdn
   </variableNames>
....
```

TPQW forcings file:
** Note, the default case is looking for specific humidity `QBOT`, but relative humidity `RH` (rh) is also fine, so long as this variable name is provided
```xml
<?xml version="1.0"?>
<file id="stream" version="1.0">
....
<fieldInfo>
   <variableNames>
      ZBOT     z
      TBOT     tbot
      WIND     wind
      RH       rh
      PSRF     pbot
      FLDS     lwdn
   </variableNames>
....
```

| <span> |
| :--- |
| **NOTE:** The latitude and longitude values of the domain file, surface datasets, and your netcdf forcing files *must* match! |
| <span> |



## 5) Get the modified parameter files

The parameter files must be modified to better simulate the mountain tundra of Niwot Ridge. The parameter changes are summarized in the table below (After Wieder et al. 2017)

| Parameter                              | name in model     | default value in clm 4.5/5.0 | changed value | 
| -------------------------------------- | ----------------- | ---------------------------- | ------------- |
| Leaf C:N ratio                         | leaf_cn           | 25/28.0269058295964          | 32            |
| Fine root: leaf ratio                  | froot_leaf        | 1/1.5                        | 2             |



If you are using parameter files that have already been modified, copy them into the case directory:

```bash
cp ../Surface_and_parameter_data/clm5_params.c200519_modfrootleaf_modleafcn.nc /path/to/your/case/case_dir
```

If you need to implement the modifications yourself follow the instructions below:

**To customize parameter files yourself:**  

1) Find location of default CLM parameter file:

```bash
grep "param" /path/to/your/case/case_dir/Buildconf/clmconf/lnd_in
> paramfile = '/glade/p/cesmdata/cseg/inputdata/lnd/clm2/paramdata/clm5_params.c200519.nc'
```
2) Copy the original parameter file to your case directory

```bash
cp /glade/p/cesmdata/cseg/inputdata/lnd/clm2/paramdata/clm5_params.c200624.nc /path/to/your/case/case_dir/clm5_params.c200624.nc
```

3) Use `nco` commands to change parameter files. In the parameter files change
```bash
module load nco # load nco module if not already loaded
# Modify the fine root:leaf ratio to 2
# Modify the leaf C:N ratio to 32, 
# the easiest way to do this is to change the value for all parameters
ncap2 -O -s 'froot_leaf=froot_leaf * 0+2;leafcn=leafcn * 0+32' clm5_params.c200624.nc clm5_params.c200624_NWTmods.nc
```

## 6) Change the user namelist files
The user namelist files are used to specify modifications to the run that are not already specified in the parameter files or *.xml files. We will use them to specify the locations of our modified surface dataset, and parameter files, as well as specifying the kind of output we want.

1) First, set the surface dataset file. This file differs for the DM, MM, WM, and SB communities. 

   - Several changes have been made to the surface datasets. All have had their clay/sand percentages adjusted to mimic the rocky soil more accurately. 
   - In addition, a different soil depth must be used for each vegetation community: 


| Community         | Soil Depth (cm) | Snow Conditions (% relative to obs) [these modifications are in the forcing netcdf files] | Surface file name |
| ----------------- | --------------- | ----------------------------------- | ----------------- | 
| Fellfield (FF)    | 70              | 10, but 25 March, April, May (MAM)  | `surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_70cm_soildepth.nc`  |
| Dry meadow (DM)   | 100             | 10, but 25 March, April, May (MAM)  |`surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_100cm_soildepth.nc` |
| Moist meadow (MM) | 130             | 100                                 | `surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_130cm_soildepth.nc` |
| Wet meadow (WM)   | 100             | 75, + runoff from simulation of MM  |`surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_100cm_soildepth.nc` |
| Snowbed (SB)      | 70              | 200                                 | `surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_70cm_soildepth.nc`  |
| | | | (Table adapted from Wieder et al. 2007) |


2) Use the following code to set the surface dataset file. Make sure you are changing the surface dataset depending on the vegetation community you are using.

```bash
echo "fsurdat= '/path/to/your/surfacedata/surfdata_0.9x1.25_16pfts_CMIP6_simyr1850_254.42_40.05_c170706_pctsandclay_100cm_soildepth.nc'" >> user_nl_clm
```

3) Set the location of the modified parameter file. This file contains important changes to leafcn and frootleaf; it must be in the case directory
```bash
echo "paramfile= '/path/to/your/case/case_dir/clm5_params.c200519_modfrootleaf_modleafcn.nc'" >> user_nl_clm
```

4) Next, set the `user_nl_datm` files so that the atmosphere cycles and uses nearest neighbor cells for deposition.
```bash
# Tell the atmosphere that when in doubt it should use the closest neighbor cell to our location for N deposition/aerosols etc.
echo "mapalgo= 'nn','nn','nn','nn','nn'" >> user_nl_datm
echo "taxmode= 'cycle','cycle','cycle','cycle','cycle'" >> user_nl_datm
```

5) Finally, set the output variables and frequency in user_nl_clm; Since this is a spinup run, we will write less frequently to the output files. And writeout fewer variables.
```bash
echo "hist_mfilt= 20" >> user_nl_clm # 20*8760 hours of data will be saved in each output file; i.e. 10 years per output file. 
echo "hist_nhtfrq= -8760" >> user_nl_clm # number of hours in a year; 

# These lines reduce the amount of times .h0. files are written.
echo "hist_empty_htapes= .true." >> user_nl_clm
# What history variables should be recorded from simulation
echo "hist_fincl1  ='TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC', 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO'" >> user_nl_clm
```

## **STOP:** Clone Case Option 
If you would like to save time by cloning your case to create other vegetation community simulations, this is the best point to do so.

See: [Cloning your case](#cloning-your-case) for details. 

## 7) Run the case in Advanced Decomposition mode
Because you have set `CLM_ACCELERATED_SPINUP` to `on`, the run will be in advanced decomposition mode. 

```bash
./case.build
./case.submit
```

## 8) Check Spinup

When the case has finished running, copy the spinup script from the CLM code base into the case directory and modify its options.

```bash
# location of script:
cp ~/clm5.0-master/tools/contrib/SpinupStability.ncl /path/to/your/case/case_dir
```
Modify script: 

```bash
;  SPT (single point) NWT BGC 2008-2012
  caseid = "case_name" ; CHANGE ME to case name
  username = "USER"   ; CHANGE ME to your user name
  annual_hist = True
  region = "SPT"                 ; Global, Arctic, or SPT (single point)
  subper = 20                    ; Subsampling period in years; CHANGE ME to # years forcing data that you are useing, (or the frequency of your history file output).
  
```
Optional: ask the script to generate spinup plots. 

```bash
 do_plot = True
```

Run the script
```bash
module load ncl # load the NCAR command language module
ncl SpinupSustainability.ncl # Run the script
```
If the spinup looks good, proceed. Otherwise, run for longer until it is spun up or troubleshoot the model if things look fishy.

## 9) Run in PostAD mode

Turn off advanced decomposition mode for 200 years and rerun. 

```bash
./xmlchange CLM_FORCE_COLDSTART=off 
./xmlchange CLM_ACCELERATED_SPINUP=off 
```

Set the simulation to run for 200 years starting from where it left off before
- turn off AD and coldstart

```bash
./xmlchange CLM_FORCE_COLDSTART=off
./xmlchange CLM_ACCELERATED_SPINUP=off

./xmlchange STOP_N=200
./xmlchange REST_N=100
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=0
./xmlchange RUN_REFDATE=0408-01-01
./xmlchange RUN_STARTDATE=0408-01-01
```

Tell the model where to find the proscribed conditions with which to start the model.
```bash
echo "finidat= '/glade/scratch/$USER/archive/case_name/rest/0408-01-01-00000/case_name.clm2.r.0408-01-01-00000.nc'" >> user_nl_clm
```

Resubmit the case, this time accellerated decomposition is off. 
```bash
./case.submit
```

## 10) Check the spinup again

Run the script

```bash
module load ncl # if not already loaded, load ncl module
ncl SpinupSustainability.ncl
```
If spinup looks good, proceed with transient run.

## 11) Run transient run with forcings

A transient run, refers to a run in which the actual experimental data is generated. Usually it takes place under experimental conditions as well. 

Point to new postAD restart file in user_nl_clm
```bash
echo "finidat= '/glade/scratch/$USER/archive/case_name/rest/0608-01-01-00000/case_name.clm2.r.0608-01-01-00000.nc'" >> user_nl_clm
```

Change model output to half-hourly data in user_nl_clm and change output variables to produce variables that are comparable to the observations. Below:

- **Older entries in user_nl_clm used for spinup (e.g., `finidat`, `hist_mfilt`, etc.) should be commented out or deleted**
- h0 = standard monthly means
- h1 files will be annual, 30 minute output

```bash
echo "finidat= '/glade/scratch/wwieder/archive/clm50bgc_NWT_ff/rest/2608-01-01-00000/clm50bgc_NWT_ff.clm2.r.2608-01-01-00000.nc' " >> user_nl_clm
echo "hist_nhtfrq = 0,1" >> user_nl_clm
echo "hist_mfilt  = 1,17520" >> user_nl_clm
echo "hist_fincl2  = 'SNOWICE','FSDS','FLDS','FSR','FSA','FIRE','FIRA','FSH','FCTR','FCEV','FGEV','FGR','FGR12','FSM','TSOI','COSZEN','RAIN','SNOW','H2OSOI','FPI','HR','ELAI','BTRAN2','FPSN','AR','RSSUN','RSSHA','FSH_G','RHAF','RH_LEAF','RH','T10','TG','SABG','SABV','EFLX_LH_TOT','SNOW_DEPTH','SOILLIQ','TV','QRUNOFF','INT_SNOW','GPP','NPP','TOTVEGC','TOTSOMC','TOTECOSYSC','WOODC','NEE','AGNPP','BGNPP','FPI'" >> user_nl_clm
```

Set env_run.xml up for transient run:

```bash
./xmlchange STOP_N=10 
./xmlchange RUN_STARTDATE=2008-01-01 
./xmlchange DATM_CLMNCEP_YR_END=2017 

```

You'll also have to make sure the **`user_datm.streams.txt.CLMGSWP3v1`** files have years for the full study period (2008-2017 in this case). I just added these manually.

Check that your namelist changes work
```bash
./preview_namelists
```

Submit the case

```bash
./case.submit
```

# Cloning your case

If you have created one case for a particular vegetation community, you can easily clone that case to save time in the setup of the other vegetation communities. To clone a case:

```bash
cd /glade/u/home/$USER/clm5.0-master/cime/scripts

./create_clone --clone /path/to/your/case/new_case_dir --case /path/to/your/case/case_dir
```

If you are cloning the case to alter vegetation community settings make sure you have changed the following before running the new case: 

 - [ ] user_datm streams must point to the correct forcing data
 - [ ] surface dataset files must be correct for the vegetation community
 - [ ] make sure your user_nl_clm is still pointing to the right surface dataset, parameter file, etc. (if you want you can copy the parameter file into the new case directory, but need to point to it in user_nl_clm)

# After you're done

The `flow.sim.R` script expects 30 minute data. Concatenating the data into a single file makes this easier.  If you followed the instructions above, the code below will concatenate your *.h1.* files together

```bash
cd /scratch/$USER/archive/$CASE_NAME/lnd/hist
module load nco
ncrcat *.h1.* CASE_NAME.clm2.h1.2008-2017.nc
```

Then, from a new terminal window, you can scp the file onto your local machine

```bash
scp $USER@cheyenne.ucar.edu:/glade/scratch/$USER/archive/$CASE_NAME/lnd/hist/$CASE_NAME.clm2.h1.2008-2017.nc .
```


### Now we can finally look at the model and observations using the `flow.obs.R` and `flow.sim.R` scripts

