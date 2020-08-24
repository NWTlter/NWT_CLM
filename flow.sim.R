# Check spinup

# Plotting basic graphs to check spinup of model
rm(list = ls())
#source("~/Documents/Fierer_lab/NCAR/NEON_LTER_NCAR/NWT_analyses/Tvan_plots/CLM_utlis.R")

##############################################################################
# Dependencies
##############################################################################

#Call the R HDF5 Library
packReq <- c("ncdf4","dplyr", "tidyr", "lubridate")

#Install and load all required packages
lapply(packReq, function(x) {
  print(x)
  if (require(x, character.only = TRUE) == FALSE) {
    install.packages(x)
    library(x, character.only = TRUE)
  }})

##############################################################################
# Workflow parameters
##############################################################################
#### Output Options ####
# Base directory for output
DirOutBase <- paste0("~/Downloads/SIM/data/")

#### Input options ####
# The input directory where simulation data is located
DirIn <- "~/Downloads/"

# The name of the netcdf file from the simulation you want to work with
ncdf_fp <- "2000datm_CLM50bgc_nwt_DM.clm2.h0.2008-01-01-00000.nc"
ncdf_fp <- "test.nc"

#### Vegetation Community ####
# Which vegetation community is this simulation for?
veg_com <- "DM" # Options: "FF", "DM", "WM", "MM", "SB", NA

##############################################################################
# Static workflow parameters - these are unlikely to change
##############################################################################
# Extract a title from netcdf name
title <- sub("\\.nc$","", basename(ncdf_fp))

# Output subdirector is the DirOutBase + the title of the netcdf
DirOut <- paste0(DirOutBase, title)

# Create output directory if it doesn't exist
if (!dir.exists(DirOut)) dir.create(DirOut, recursive = TRUE)



# dir <- "~/Downloads/"
# title_list <- c("2000datm_CLM50bgc_nwt.clm2.h0.200yr_spinup",
#                 "2000datm_CLM50bgc_nwt.clm2.h0.100yr_spinup_frootleaf_parameter_changes",
#                 "2000datm_CLM50bgc_nwt.clm2.h0.150yr_spinup_modleafcn_parameter_changes",
#                 "2000datm_CLM50bgc_nwt.clm2.h0.400yr_spinup_leafcn_frootleaf_sandclay_changes",
#                 "2000datm_CLM50bgc_nwt_DM.clm2.h0.2008-01-01-00000",
#                 "2000datm_CLM50bgc_nwt_newTvan_base.clm2.h0.monthly_spinup", 
#                 "2000datm_CLM50bgc_nwt_newTvan_DM.clm2.h0.monthly_100yrs",
#                 "2000datm_CLM50bgc_nwt_newTvan_base.clm2.h0.monthly_1850_50yrs_test")

# Choose simulation variation
# 1 = code changes, nothing else - run for 200 yrs
# 2 = code changes plus frootleaf changed from 1.5 to 2 - run for 100 yrs
# 3 = code changes plus leafcn changed from 28.026ish to 32 - run for 150 yrs
# 4 = code changes plus frootleaf, leafcn, and 39%sand23%clay soil type
# 5 = code changes plus frootleaf, leafcn, and 39%sand23%clay soil type and dry meadow surface file; also a fully spun-up transient run.
# 6 = monthly measurements from tvan 2008-2016 forcings; ad mode code changes plus frootleaf, leafcn, and 39%sand23%clay soil type
# 7 = monthly measurements from tvan 2008-2016 forcings; ad mode code changes plus frootleaf, leafcn, and 39%sand23%clay soil type and dry meadow surface file; 100 yrs ad mod
# title <- title_list[5]
################################################################################
# Helper functions - for downloading and loading data
################################################################################
#---------------read in CLM variables----------------------------------
extract_CLM_vars <- function(infile, nsteps, vars) {
  # Extract variables from netcdf history files
  # file = filepath to netcdf file
  # nsteps = how many 30minute timesteps worth of data are there?
  # vars = the variables to extract from the netcdf file
  require(ncdf4)
  
  
  #infile   <- "~/Downloads/test.nc"
  
  Data.clm <- nc_open(infile)   
  nsteps   <- 60 #48 * (365*5 + 2) #48 * (365*6 + 2)
  print(paste("The file has",Data.clm$nvars,"variables"))
  print(paste("The variables are:"))
  print(paste(names(Data.clm$var)))
  summary(Data.clm)
  #print(Data.clm)
  
  # vars <- c("mcdate", "mcsec", "DZSOI", "TV", "T10", "RH", "RH_LEAF", "BTRAN", "FSH",
  #           "EFLX_LH_TOT", "FGR", "FGR12", "FSDS", "FIRA", "FPSN", "GPP", "NEE", "ELAI",
  #           "FSR", "FSA", "TSOI", "SOILLIQ", "H2OSOI", "SNOW", "RAIN", "SNOW_DEPTH",
  #           "QVEGT", "QVEGE", "QSOIL", "FCEV", "FCTR", "FGEV")
  
  
  # If vars is not given, set it to be all the vars in the file
  # otherwise, take the user-supplied vars
  if(missing(vars)) {
    vars <- names(Data.clm$var)
  } else {
    # Check if vars exist in output, and throw a warning if they don't
    if (any(!(vars %in% names(Data.clm$var)))) {
      warning(paste0("Cannot load all variables because some do not exist. \n",
                     "The following variables do not exist in your history file: \n",
                     paste(vars[which(!(vars %in% names(Data.clm$var)))], collapse = ", ")))
    }
  }
  
  
  
  # list of variables pulled from netcdf
  varlist <- vector(length = length(vars), "list")
  names(varlist) <- vars
  # list of dimensions of those variables
  vardim <- vector(length = length(vars), "list")
  names(vardim) <- vars
  # list of the units for each variable
  varunitlist <- vector(length = length(vars), "list") 
  names(varunitlist) <- vars
  # list of the long name of each variable
  varlongname <- vector(length = length(vars), "list")
  names(varlongname) <- vars
  
  # Load variables and their units into lists; Separate out soil layer and other
  # layered variables into separate columns
  if (exists("varord")) {rm(varord)}
  writeLines("Loading variables into list...")
  for (i in seq_along(vars)) {
    print(vars[i])
    try(vardim[[i]] <- dim(ncvar_get(Data.clm, vars[i])))
    # if the variable has 2 dimensions, flatten them into separate columns
    if (length(vardim[[i]]) == 2 ) {
      print(paste0("Variable ", vars[i], " has ", vardim[[i]][1], " columns."))
      try(tmp_var <- ncvar_get(Data.clm, vars[i]))
      try(tmp_varunit <- ncatt_get(Data.clm, vars[i], "units"))
      try(tmp_varlongname <- ncatt_get(Data.clm, vars[i], "long_name"))
      for (j in 1:vardim[[i]][1]) {
        tmp_name <- paste0(vars[i], "_", j)
        tmp_longname <- tmp_varlongname
        varord <- c(varord, tmp_name)
        try(varlist[[tmp_name]] <- tmp_var[j,])
        try(varunitlist[[tmp_name]] <- tmp_varunit)
        if(tmp_longname$hasatt) {
          tmp_longname$value <- paste0(tmp_longname$value, 
                                       " at level ", j)
          # print(paste0("tmp_varlongname has attribute, it is ", tmp_longname$value))
        } else {
          tmp_longname$value <- NA
        }
        try(varlongname[[tmp_name]] <- tmp_longname)
      }
    } else if (exists("varord")) { # otherwise simply add the dimension to the list
      #print(paste0("this is the regular loop for 1 dimensional variables"))
      varord <- c(varord, vars[i])
      try(varlist[[i]] <- ncvar_get(Data.clm, vars[i]))
      try(varunitlist[[i]] <- ncatt_get(Data.clm, vars[i], "units"))
      try(varlongname[[i]] <- ncatt_get(Data.clm, vars[i], "long_name"))
    } else { # for the first variable create the varord vector
      #print("this should only happen once, for the first variable")
      varord <- vars[i]
      try(varlist[[i]] <- ncvar_get(Data.clm, vars[i]))
      try(varunitlist[[i]] <- ncatt_get(Data.clm, vars[i], "units"))
      try(varlongname[[i]] <- ncatt_get(Data.clm, vars[i], "long_name"))
    }
  }
  
  # reorder variables
  varlist <- varlist[varord]
  varunitlist <- varunitlist[varord]
  varlongname <- varlongname[varord]
  
  # Optional creation of derived variables if they exist:
  # Create net radiation if FSA and FIRA are present
  if (("FSA" %in% varord) & ("FIRA" %in% varord)) {
    varlist$RNET <- varlist$FSA - varlist$FIRA
    varunitlist$RNET <- varunitlist$FSA
    varlongname$RNET <- list(hasatt = TRUE, value = "net radiation: FSA-FIRA")
  }
  
  nc_close(Data.clm)
  # print(mean(GPP)*3600*24*365)
  # print(mean(LH))
  # case[e]
  return(list(data = varlist, units = varunitlist, longname = varlongname))
  
  
}

summarize_vars_by_time <- function(var, unitlist, ncdata, veg_com = NA) {
  library(dplyr)
  library(lubridate)
  library(tidyr)
  if ( missing(unitlist) ) {
    unitlist <- rep("", length(var))
    names(unitlist) <- var
  }
  
  var.plot.all <- as.data.frame(ncdata[c("mcsec", "mcdate", var)]) %>%
    mutate(mcdate = stringr::str_pad(mcdate, width = 8, side = "left", pad = "0"),
           date_UTC = as.Date(mcdate, format = "%Y%m%d"), 
           timestamp_UTC = as.Date(date_UTC, origin = paste0(date, " 00:00:00")) +
             lubridate::seconds(mcsec)) %>%
    # Convert the timezone
    mutate(timestamp = lubridate::with_tz(timestamp_UTC, tzone = "MST"),
           date = as.Date(timestamp, tz = "MST"),
           Hour = lubridate::hour(timestamp) + lubridate::minute(timestamp)/60,
           month = as.factor(lubridate::month(date)),
           year = as.factor(lubridate::year(date)),
           DoY = lubridate::yday(date)) %>%
    select(mcsec, mcdate, timestamp, date, timestamp_UTC, date_UTC, year, month,
           DoY, Hour, all_of(var)) %>%
    mutate(ObsSim = "Sim") %>%
    mutate(veg_com = veg_com)
  
  var.plot.diurnal <- var.plot.all %>%
    mutate(MonGroup = ifelse(month %in% c(12,1,2), "DJF",
                             ifelse(month %in% c(3,4,5), "MAM", 
                                    ifelse(month %in% c(6,7,8), "JJA", "SON")))) %>%
    group_by(MonGroup, Hour) %>% 
    mutate(across(all_of(var), .fns = list(houravg = ~mean(., na.rm = TRUE),
                                           hoursd  = ~sd(., na.rm = TRUE)))) %>%
    ungroup() %>% 
    select(MonGroup, Hour, contains("houravg"), contains("hoursd")) %>%
    unique() %>%
    mutate(ObsSim = "Sim") %>%
    mutate(veg_com = veg_com)
  
  # var.plot.diurnal %>%
  #   select(MonGroup, Hour, ObsSim, starts_with(c("RNET", "FSH", "EFLX", "GPP"))) %>%
  #   pivot_longer(ends_with("houravg"), names_to = "MeanMetric",
  #                values_to = "MeanValue") %>%
  #   pivot_longer(ends_with("hoursd"), names_to = "SDMetric",
  #                values_to = "SDValue") %>%
  #   mutate(MeanMetric = gsub("_houravg", "", MeanMetric),
  #          SDMetric = gsub("_hoursd", "", SDMetric)) %>%
  #   filter(MeanMetric == SDMetric) %>%
  #   ggplot(aes(x = Hour)) +
  #   geom_ribbon(aes(ymin = MeanValue-SDValue, ymax = MeanValue+SDValue),
  #               alpha = 0.5) +
  #   geom_line(aes(y = MeanValue, color = ObsSim)) +
  #   facet_wrap(MonGroup~MeanMetric)
  
  
  var.plot.doy <- var.plot.all %>%
    group_by(DoY) %>% 
    mutate(across(all_of(var), .fns = list(doyavg = ~mean(., na.rm = TRUE),
                                           doysd  = ~sd(., na.rm = TRUE)))) %>%
    ungroup() %>%
    select(DoY, contains("doyavg"), contains("doysd")) %>%
    unique() %>%
    mutate(ObsSim = "Sim") %>%
    mutate(veg_com = veg_com)
  
  var.plot.ann <- var.plot.all %>%
    group_by(year) %>% 
    mutate(across(all_of(var), .fns = list(yearavg = ~mean(., na.rm = TRUE),
                                           yearsd  = ~sd(., na.rm = TRUE)))) %>%
    ungroup() %>%
    select(year, contains("yearavg"), contains("yearsd")) %>%
    unique() %>%
    mutate(ObsSim = "Sim") %>%
    mutate(veg_com = veg_com)
  
  ##### ####  
  #   # get an average of every month
  #   group_by(month, year) %>%
  #   mutate(across(all_of(var), .fns = list(monavg = ~mean(., na.rm = TRUE)))) %>%
  #   ungroup() %>%
  #   # get a daily average
  #   group_by(DoY, month, year) %>%
  #   mutate(across(all_of(var), .fns = list(dailyavg = ~mean(., na.rm = TRUE)))) %>%
  #   ungroup()
  # 
  # var.plot.timestep <- var.plot.all %>%
  #   select(!contains("monavg"), !contains("dailyavg")) %>%
  #   pivot_longer(cols = all_of(var), names_to = "Variables", values_to = "Values") %>%
  #   mutate(Variables_Units = unitlist[Variables])
  # 
  # # Monthly averages
  # var.plot.month <- var.plot.all %>%
  #   select(mcsec, mcdate, date, month, year, DoY, Hour, timestamp, timestamp_UTC, 
  #          contains("monavg")) %>%
  #   filter(grepl(".{1,4}-.{2}-01", timestamp_UTC)) %>% 
  #   filter(mcsec == 0) %>%
  #   pivot_longer(cols = contains("monavg"), 
  #                names_to = "Variables", values_to = "Values") %>%
  #   mutate(Variables = gsub("_monavg", "", Variables)) %>%
  #   mutate(Variables_Units = unitlist[Variables])
  # 
  # # Daily averages
  # var.plot.daily <- var.plot.all %>%
  #   select(mcsec, mcdate, date, month, year, DoY, Hour, timestamp, 
  #          contains("dailyavg")) %>%
  #   filter(mcsec == 0) %>%
  #   pivot_longer(cols = contains("dailyavg"), 
  #                names_to = "Variables", values_to = "Values") %>%
  #   mutate(Variables = gsub("_dailyavg", "", Variables)) %>%
  #   mutate(Variables_Units = unitlist[Variables])
  #### #### 
  return(list(diurnal_seasonal = var.plot.diurnal, 
              daily = var.plot.doy,
              annual = var.plot.ann,
              all_wide = var.plot.all))
}
################################################################################
# Extract the CLM variables into an R-friendly format
################################################################################
nc <- extract_CLM_vars(infile = paste0(DirIn, title, ".nc"))


# Subset nc data to only include the time step data
# get list of variables that are not recorded each time step
run_data_names <-  sapply(nc$data, function(i) length(i) < length(nc$data$nstep))
run_data <- nc$data[!run_data_names]

# list of variable units
run_units <-  sapply(names(run_data), function(x) ifelse(nc$units[[x]]$hasatt,
                                                         nc$units[[x]]$value,
                                                         NA))
run_units <- data.frame(t(run_units), stringsAsFactors = FALSE)

# list of variable descriptions
run_longname <- sapply(names(run_data), function(x) nc$longname[[x]]$value)
run_longname <- data.frame(t(run_longname), stringsAsFactors = FALSE)

################################################################################
# Unit conversions
################################################################################
# Convert all Kelvin temperatures to Celsius
# find variables with "temperature" in the longname, then subset those to variables
# with units of K
temperature_vars <- grep("temperature", run_longname) 
temperature_vars <- temperature_vars[grepl("K", run_units[temperature_vars])]

writeLines("Converting Kelvin variables to Celsius")
for (i in seq_along(temperature_vars)) {
  #print(run_longname[temperature_vars[i]])
  #print(run_units[temperature_vars[i]])
  run_data[[temperature_vars[i]]] <- run_data[[temperature_vars[i]]] - 273.15
  run_units[temperature_vars[i]] <- "C"
}

# Convert all snow depth into cm
writeLines("Converting SNOW_DEPTH from m to cm")
run_data[["SNOW_DEPTH"]] <- run_data[["SNOW_DEPTH"]] * 100 # from m to cm
run_units["SNOW_DEPTH"] <- "cm"

# Convert soil moisture into percent
soilmoist_vars <- grep("H2OSOI", names(run_longname))

writeLines("Converting soil moisture from ratio to %")
for (i in seq_along(soilmoist_vars)) {
  run_data[[soilmoist_vars[i]]] <- run_data[[soilmoist_vars[i]]] * 100
  run_units[soilmoist_vars[i]] <- "% mm/mm"
}

################################################################################
# Summarize clm variables at different time chunks (i.e. diurnally, hourly etc.)
################################################################################
writeLines(paste0("Selecting variables to compare to obs and creating diurnal, ",
                  "daily, and annual summaries of them"))

# Select variables to pull from model
var <- c("FSH", # sensible heat flux (W/m^2)
         "T10", # temperature at 2m (C)
         "RNET", # net radiation FSA-FIRA (W/m^2)
         "TSOI_2", # Soil temperature (2nd layer; 4cm, 2-6cm layer) (C)
         "TSOI_3", # Soil temperature (3rd layer; 9cm, 6-12cm layer) (C)
         "TSOI_5", # Soil temperature (5th layer; 26cm, 20-32cm layer) (C)
         "H2OSOI_2", # volumetric soil moisture (2nd layer; 4cm, 2-6cm layer) (mm/mm)
         "H2OSOI_3", # volumetric soil moisture (3rd layer; 9cm, 6-12cm layer) (mm/mm),
         "H2OSOI_5", # volumetric soil moisture (5th layer; 26cm, 20-32cm layer) (mm/mm)
         "SNOW_DEPTH", # snow height of snow covered area (cm)
         "EFLX_LH_TOT", # total latent heat flux [+ to atm], (W/m^2)
         "GPP", # Gross primary production (gC/m^2/s)
         "AGNPP", # Aboveground net primary productivity (gC/m^2/s)
         "NPP") # Net primary production(gC/m^2/s))

# Select the appropriate soil layers for comparison
if (veg_com == "FF") {
  var <- var[!grepl("SOI.{1,}2", var)] # FF is measured at 10cm depth
} else {
  var <- var[!grepl("SOI.{1,}3", var)] # all other communities measured at 5cm
}



# get the units associated with the variables
unitlist <- run_units[var]

writeLines("Selected Variables are: ")
for (i in seq_along(var)) {
  variable <- var[i]
  writeLines(paste0(variable, ": ", run_longname[variable], " (",
                    run_units[variable], ")"))
}

# Extract the three time summaries of the variables
nc_format <- summarize_vars_by_time(var = var,
                       unitlist = unitlist, 
                       ncdata = run_data,
                       veg_com = veg_com)



################################################################################
# Write out Simulation data
################################################################################
writeLines("Writing out simulation data")

# Write out diurnal summaries for each season
write.table(nc_format$diurnal_seasonal, 
            file = paste0(DirOut, "/Diurnal_seasonal_summaries_", title, ".txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Write out daily summaries for each DoY
write.table(nc_format$daily, 
            file = paste0(DirOut, "/Daily_summaries_", title, ".txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Write out annual summaries for each year
write.table(nc_format$daily, 
            file = paste0(DirOut, "/Annual_summaries_", title, ".txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Write out 30-minute data 
write.table(nc_format$all_wide, 
            file = paste0(DirOut, "/Unsummarized_30min_", title, ".txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Write out unit definitions
write.table(unitlist, 
            file = paste0(DirOut, "/Unit_Definitions_", title, ".txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")


# #####################################################################################
# diag.plots <- plot_var(var = var, unitlist = unitlist,
#          ncdata = run_data, strtyr = 0, endyr = 20)
# 
# pdf(paste0(dir, title, "_diag_plots.pdf"))
# diag.plots
# dev.off()
# 
# # Read in flux tower data:
# library(REddyProc)
# # data.G <- fLoadTXTIntoDataframe("~/Downloads/tvan_obs/Tvan_flux_OBS_2008-2013b.txt")
# # data.G
# nsteps <- 48 * (365*6 + 2)
# 
# data.sno <- read.csv("~/Downloads/tvan_obs/SnowDepth_daily.csv")
# names(data.sno)
# 
# Data.flx <- fLoadTXTIntoDataframe("~/Downloads/tvan_obs/Tvan_flux_OBS_2008-2013b.txt")
# names(Data.flx)
# #  units(Data.flx)
# maxYear <- max(Data.flx$Year)
# 
# #read in ground flux measurements:
# dir2 <- '/Users/wwieder/Desktop/Working_files/Niwot/NR_fluxes/TVan/G/'
# Data.G  <- fLoadTXTIntoDataframe("~/Downloads/tvan_obs/Wieder_Niwot_CLM_G.txt")
# Gobs    <- Data.G$G[Data.G$Year <= maxYear]
# 
# 
# # Prepare data for comparison plots - need, timestamp + variables of interest
# # plus any units need to be converted to be the same
# # Flux tower data
# Data.flx.plot <- Data.flx %>%
#   mutate(timestamp = as.POSIXct(paste0(Year, "-", MO, "-", DD, " ", HR, ":", MM), 
#                                 tz = "MST")) %>%
#   left_join(data.sno %>%
#               mutate(timestamp = as.POSIXct(paste0(year_mean, "-", 
#                                                    mo_mean, "-", 
#                                                    day_mean, 
#                                                    " ", "00:00"), 
#                                             tz = "MST")),
#             by = c("Year" = "year_mean", "MO" = "mo_mean", "DD" = "day_mean")) %>%
#   rename(date = timestamp.y, timestamp = timestamp.x) %>%
#   select(Year, MO, DD, HR, MM, DecimalDate, timestamp, date, everything()) %>%
#   group_by(Year, MO) %>%
#   mutate_at(vars(NEE:SnoDep_med), ~mean(., na.rm = TRUE)) %>%
#   ungroup() %>%
#   filter(DD == "1") %>%
#   select(-HR, -MM, -DecimalDate, -timestamp, -date) %>%
#   unique() %>%
#   mutate(timestamp = as.POSIXct(paste0(Year, "-",
#                                        MO, "-",
#                                        DD, " ", "00:00"),
#                                 tz = "MST")) %>%
#   select(-Year, -MO, -DD) %>%
#   rename(GPP = GPP_f, Rn = Rn_MDS) %>%
#   select(timestamp, everything())
# 
# run_data.plot <- as.data.frame(run_data) %>%
#   select(mcdate, mcsec, AGNPP:WOODC) %>%
#   mutate(mcdate = stringr::str_pad(mcdate, width = 8, side = "left", pad = "0"),
#          time = as.POSIXct(lubridate::ymd_hm(paste0(mcdate, "0000"), tz = "UTC") + 
#            lubridate::dseconds(mcsec), tz = "UTC"),
#          timestamp = as.POSIXct(lubridate::with_tz(time,tzone = "MST"))
#          ) %>%
#          #timestamp = as.POSIXct(paste0(mcdate," 00:00"), 
#          #                        format = "%Y%m%d HH:MM", tz = "GMT")) %>%
#   select(timestamp, AGNPP:WOODC) %>%
#   # modify some key variables
#   mutate(Rn = FSA - FIRA,
#          T10 = T10 - 273.15, # convert to deg C
#          GPP = GPP * 1000, # convert g to mg
#          NEE = NEE * 1000, # convert g to mg
#          SNOW_DEPTH = SNOW_DEPTH * 100) # convert m to cm
# 
# 
# sim_var <- c("NEE", "GPP", "EFLX_LH_TOT", "FSH", "T10", "Rn", "SNOW_DEPTH", "SOILLIQ_3")
# obs_var <- c("NEE", "GPP", "LE", "H", "Tair", "Rn", "SnoDep_med", "SoilMoisture")
# sim_data <- run_data.plot
# obs_data <- Data.flx.plot
# unitlist <- c("mgC/m^2/s", "mgC/m^2/s", "W/m^2" , "W/m^2", "deg C", "W/m^2", "cm", "kg/m2")
# names(unitlist) <- obs_var
# 
# 
# obs_comp_plot <- plot_var_comp(sim_var = sim_var, obs_var = obs_var, sim_data = sim_data, obs_data = obs_data, unitlist = unitlist, title = title)
# 
# pdf(paste0(dir, title, "_obs_comp.pdf"))
# obs_comp_plot
# dev.off()
