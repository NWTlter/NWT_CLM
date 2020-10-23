rm(list = ls())

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
# Base directory for output, to conform with plotting script
user = 'wwieder'
if (user ==  'wwieder') {
  DirBase <- "~/Desktop/Working_files/Niwot/CLM/"
  DirOutBase <- paste0(DirBase,"SIM/")
} else {
  DirBase <- "~/Downloads/"
  DirOutBase <- paste0(DirBase,"SIM/")
}

# Case name of the simulation to create an output subdirectory (optional), if you don't want
# to specify a case name, set equal to ""
# This is useful if you are running many cases
case_name <- "clm50bgc_NWT"

#### Input options ####
# The input directory where simulation data is located
DirIn <- paste0(DirOutBase,"clm_history_files/")

# The names of the netcdf files from each simulation you want to work with.
# If you don't have a netcdf for a particular vegetation community, leave
# the file path blank.
ff_ncdf_fp <- "clm50bgc_NWT_ff_lowSAT.clm2.h1.2008-2017.nc" # fell field
dm_ncdf_fp <- "clm50bgc_NWT_dm_lowSAT.clm2.h1.2008-2017.nc" # dry meadow
wm_ncdf_fp <- "" # wet meadow
mm_ncdf_fp <- "clm50bgc_NWT_mm.clm2.h1.2008-2017.nc" # moist meadow
sb_ncdf_fp <- "clm50bgc_NWT_sb.clm2.h1.2008-2017.nc" # snow bed

#### Extra Variable Choice ####
# The names of any optional extra variables the user would like to extract
# WARNING: No conversions are made to the units of user-specified optional 
# variables
# 
# Variables that are extracted by default are:
#      "FSH", # sensible heat flux (W/m^2)
#      "T10", # temperature at 2m (C)
#      "RNET", # net radiation FSA-FIRA (W/m^2)
#      "TSOI_2", # Soil temperature (2nd layer; 4cm, 2-6cm layer) (C)
#      "TSOI_3", # Soil temperature (3rd layer; 9cm, 6-12cm layer) (C)
#      "TSOI_5", # Soil temperature (5th layer; 26cm, 20-32cm layer) (C)
#      "H2OSOI_2", # volumetric soil moisture (2nd layer; 4cm, 2-6cm layer) (mm/mm)
#      "H2OSOI_3", # volumetric soil moisture (3rd layer; 9cm, 6-12cm layer) (mm/mm),
#      "H2OSOI_5", # volumetric soil moisture (5th layer; 26cm, 20-32cm layer) (mm/mm)
#      "SNOW_DEPTH", # snow height of snow covered area (cm)
#      "EFLX_LH_TOT", # total latent heat flux [+ to atm], (W/m^2)
#      "GPP", # Gross primary production (gC/m^2/s)
#      "AGNPP", # Aboveground net primary productivity (gC/m^2/s)
#      "NPP" # Net primary production(gC/m^2/s))
usr_var <- c("ELAI", "TOTVEGC")

##############################################################################
# Static workflow parameters - these are unlikely to change
##############################################################################
# Create a list of netcdfs for each vegetation community
veg_coms_names <- c("FF", "DM", "WM", "MM", "SB")

ncdf_fp_list <- list(ff_ncdf_fp, dm_ncdf_fp, wm_ncdf_fp, mm_ncdf_fp, sb_ncdf_fp)

names(ncdf_fp_list) <- veg_coms_names

# Remove members of list without simulations
ncdf_fp_list <- ncdf_fp_list[ncdf_fp_list != ""]

# Output subdirector is the DirOutBase + the case_name
DirOut <- paste0(DirOutBase, case_name)

# Create output directory if it doesn't exist
if (!dir.exists(DirOut)) dir.create(DirOut, recursive = TRUE)

# Variables to summarize with mean annual summaries
mean_ann_sum_vars <- c("GPP", "NPP", "ET", "TOTVEGC")


################################################################################
# Helper functions - for downloading and loading data
################################################################################
#---------------read in CLM variables----------------------------------
extract_CLM_vars <- function(infile, vars) {
  # Extract variables from netcdf history files
  # file = filepath to netcdf file
  # vars = the variables to extract from the netcdf file
  require(ncdf4)

  Data.clm <- nc_open(infile) 
  # should this not be hard coded?  I'm not sure what it does later on, but it's burried in the code now?
  
  print(paste("The file has",Data.clm$nvars,"variables"))
  print(paste("The variables are:"))
  print(paste(names(Data.clm$var)))
  summary(Data.clm)

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
    # Move time back 30 minutes since timestamp seems to be off
    mutate(Hour = Hour - 0.5,
           Hour = ifelse(Hour == -0.5, 23.5, Hour)) %>%
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
  
  return(list(diurnal_seasonal = var.plot.diurnal, 
              daily = var.plot.doy,
              annual = var.plot.ann,
              all_wide = var.plot.all))
}

################################################################################
# Extract the CLM variables into an R-friendly format
################################################################################
nc_info_list <- lapply(ncdf_fp_list, function(x) {
  
  writeLines(paste0("Extracting data from ", x))
  nc <- extract_CLM_vars(infile = paste0(DirIn, x))
  
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
  
  nc_info <- list(run_data = run_data, 
                  run_units = run_units, 
                  run_longname = run_longname)
  return(nc_info)
  }
  )

################################################################################
# Unit conversions
################################################################################
nc_info_list <- lapply(nc_info_list, function(x) {
  # Convert all Kelvin temperatures to Celsius
  # find variables with "temperature" in the longname, then subset those to variables
  # with units of K
  temperature_vars <- grep("temperature", x$run_longname) 
  temperature_vars <- temperature_vars[grepl("K", x$run_units[temperature_vars])]
  
  writeLines("Converting Kelvin variables to Celsius")
  for (i in seq_along(temperature_vars)) {
    #print(x$run_longname[temperature_vars[i]])
    #print(x$run_units[temperature_vars[i]])
    x$run_data[[temperature_vars[i]]] <- x$run_data[[temperature_vars[i]]] - 273.15
    x$run_units[temperature_vars[i]] <- "C"
  }
  
  # Convert all snow depth into cm
  writeLines("Converting SNOW_DEPTH from m to cm")
  x$run_data[["SNOW_DEPTH"]] <- x$run_data[["SNOW_DEPTH"]] * 100 # from m to cm
  x$run_units["SNOW_DEPTH"] <- "cm"
  
  # Convert soil moisture into percent
  soilmoist_vars <- grep("H2OSOI", names(x$run_longname))
  
  writeLines("Converting soil moisture from ratio to %")
  for (i in seq_along(soilmoist_vars)) {
    x$run_data[[soilmoist_vars[i]]] <- x$run_data[[soilmoist_vars[i]]] * 100
    x$run_units[soilmoist_vars[i]] <- "% mm/mm"
  }
  return(x)
  })

################################################################################
# Summarize clm variables at different time chunks (i.e. diurnally, hourly etc.)
################################################################################
nc_format_list <- lapply(names(nc_info_list), function(x, opt_var = usr_var) {
  writeLines(paste0("For each vegetation community we are selecting variables to compare to obs",
                    "and creating diurnal, daily, and annual summaries of them"))
  # Set up veg_com variable based on list name
  veg_com <- x
  
  # Set up nc_info object
  nc_info <- nc_info_list[[x]]
  
  
  # Select variables to pull from model if missing above
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
  
  # add any user-specified optional variables that are not already in var
  var <- c(var, opt_var[!(opt_var %in% var)])
  
  
  # Select the appropriate soil layers for comparison
  if (veg_com == "FF") {
    var <- var[!grepl("SOI.{1,}2", var)] # FF is measured at 10cm depth (drop unused variables)
  } else {
    var <- var[!grepl("SOI.{1,}3", var)] # all other communities measured at 5cm (drop unused variables)
  }
  
  
  # get the units associated with the variables
  unitlist <- nc_info$run_units[var]

  writeLines("Selected Variables are: ")
  for (i in seq_along(var)) {
    variable <- var[i]
    writeLines(paste0(variable, ": ", nc_info$run_longname[variable], " (",
                      nc_info$run_units[variable], ")"))
  }

  # Extract the three time summaries of the variables
  writeLines(paste0("Extracting the three time summaries of the ", veg_com, " data"))
  nc_format <- summarize_vars_by_time(var = var,
                                      unitlist = unitlist,
                                      ncdata = nc_info$run_data,
                                      veg_com = veg_com)
  
  # Change the names of the soil variables to "upper" and "lower" so that FF data can
  # be combined with all other veg communities later
  
  nc_format <- lapply(nc_format, function(x) {
    names(x) <- sub("_5_", "_lower_", names(x))
    names(x) <- sub("_3_|_2_", "_upper_", names(x))
    return(x)
  })
  
  # Add the unitlist to nc_format
  nc_format$unilist <- unitlist
  
  return(nc_format)
  
})

names(nc_format_list) <- names(nc_info_list)


################################################################################
# Concatenate Vegetation Community data
################################################################################
# Create list to hold the concatenated data
nc_format <- vector(mode = "list", length = 4)
names(nc_format) <- c("diurnal_seasonal", "daily", "annual", "all_wide")

# Concatenate data
nc_format$diurnal_seasonal <- do.call(bind_rows, lapply(nc_format_list, function(x) {x$diurnal_seasonal}))
nc_format$daily <- do.call(bind_rows, lapply(nc_format_list, function(x) {x$daily}))
nc_format$annual <- do.call(bind_rows, lapply(nc_format_list, function(x) {x$annual}))
nc_format$all_wide <- do.call(bind_rows, lapply(nc_format_list, function(x) {x$all_wide}))

################################################################################
# Produce basic summaries of simulation data
################################################################################
# Mean annual summaries
mean_ann_sum <- nc_format$annual %>%
  group_by(veg_com) %>%
  select(veg_com, starts_with(mean_ann_sum_vars)) %>%
  summarize_at(vars(ends_with("_yearavg")), mean, na.rm = TRUE)

write.table(mean_ann_sum, 
            file = paste0(DirOut, "/Mean_annual_summaries_", paste(mean_ann_sum_vars, collapse = "_"), ".txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

# maximum ELAI summary
max_elai <- nc_format$all_wide %>%
  group_by(veg_com, year) %>%
  select(veg_com, year, starts_with("ELAI")) %>%
  summarize_at(vars(starts_with("ELAI")), max, na.rm = TRUE) %>%
  rename(max_ELAI = ELAI) %>%
  ungroup() %>%
  group_by(veg_com) %>%
  mutate(mean_max_ELAI = mean(max_ELAI, na.rm = TRUE))

write.table(max_elai, 
            file = paste0(DirOut, "/Max_elai_summary.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

################################################################################
# Write out Simulation data
################################################################################
writeLines("Writing out simulation data")

# Write out diurnal summaries for each season
write.table(nc_format$diurnal_seasonal, 
            file = paste0(DirOut, "/Diurnal_seasonal_summaries.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Write out daily summaries for each DoY
write.table(nc_format$daily, 
            file = paste0(DirOut, "/Daily_summaries.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Write out annual summaries for each year
write.table(nc_format$annual, 
            file = paste0(DirOut, "/Annual_summaries.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Write out 30-minute data 
write.table(nc_format$all_wide, 
            file = paste0(DirOut, "/Unsummarized_30min.txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t")

# Write out unit definitions for each vegetation community (they should all have the same units)
lapply(veg_coms_names, function(x) {
  write.table(nc_format_list[[x]]$unitlist, 
              file = paste0(DirOut, "/Unit_Definitions_", x, ".txt"),
              col.names = TRUE, row.names = FALSE, sep = "\t")
} )

# #####################################################################################
# model data have been processed for plotting & comparison with observations
# Move onto the `Obs_sim_comp_plots.R` script
# #####################################################################################

print('script complete')

