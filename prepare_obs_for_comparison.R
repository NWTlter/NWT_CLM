##############################################################################################
#' title Workflow to NCAR CLM Observational data

#' author
#' Hannah Holland-Moritz (hhollandmoritz AT gmail.com)
#' 
#' description 
#' Workflow for collating NIWOT LTER data in preparation to compare observations
#' to CLM output.


##############################################################################

##############################################################################
# Dependencies
##############################################################################
rm(list = ls())

#Call the R HDF5 Library
packReq <- c("magrittr","EML", "dplyr", "ggplot2", 
             "purrr", "tidyr", "lubridate","RCurl")

#Install and load all required packages
lapply(packReq, function(x) {
  print(x)
  if (require(x, character.only = TRUE) == FALSE) {
    install.packages(x)
    library(x, character.only = TRUE)
  }})

#Setup Environment
options(stringsAsFactors = F)

##############################################################################
#Workflow parameters
##############################################################################
#### Output Options ####
# 1) Base directory for output
# 2) Directory to download observation data to 
# 3) Location of the tvan data that was used to create forcing files
#      location of tvan data with soil information; Note: Tvan soil temperature data
#      probes from East tower do not work, so please give west tower tvan data location
# I'm trying to make it so we don't have to keep changing this...

user = 'jayka'
if (user ==  'wwieder') {
  DirOutBase <- paste0("~/Desktop/Working_files/Niwot/CLM/OBS/data")
  DirDnld = "~/Desktop/Working_files/Niwot/CLM/OBS/NWT_lter_obs_downloads"
  tvan_data_fp <- "~/Desktop/Working_files/Niwot/CLM/datav20200824T1008/data/tvan_forcing_data_precip_mods_both_towers_2007-05-11_2020-08-11.txt"
  tvan_data_soil <- "~/Desktop/Working_files/Niwot/Tvan_out_new/filtered_data/tvan_West_2007-05-09_19-00-00_to_2020-08-11_00-30-00_flux_P.csv"
  
} 
if (user ==  'jayka') {
  DirOutBase <- paste0("~/Desktop/Working_files/Niwot/CLM/OBS/data")
  DirDnld = "~/Desktop/Working_files/Niwot/CLM/OBS/NWT_lter_obs_downloads"
  tvan_data_fp <- "~/Desktop/Working_files/Niwot/CLM/data/data/tvan_forcing_data_precip_mods_both_towers_2007-05-11_2021-11-08.txt"
  tvan_data_soil <- "~/Desktop/Working_files/Niwot/Tvan_out_new/AmeriFlux_readyData/tvan_West_2007-05-09_19-00-00_to_2021-11-08_00-30-00_flux_P.csv"
  
} else { 
  DirOutBase <- paste0("~/Downloads/OBS/data") 
  DirDnld = "~/Downloads/CLM/OBS/NWT_lter_obs_downloads"
  tvan_data_fp <- "~/Downloads/CLM/datav20200816T1808/data/tvan_forcing_data_precip_mods_both_towers_2007-05-11_2020-08-11.txt"
  tvan_data_soil <- "~/Downloads/Tvan_out_new/filtered_data/tvan_West_2007-05-09_19-00-00_to_2020-08-11_07-30-00_flux_P.csv"
}

# Should a newer version of EDI data be downloaded if one is available?
getNewData = TRUE



##############################################################################
# Static workflow parameters - these are unlikely to change
##############################################################################

#Check if directory exists and create if not
if(!dir.exists(DirOutBase)) dir.create(DirOutBase, recursive = TRUE)
if(!dir.exists(DirDnld)) dir.create(DirDnld, recursive = TRUE)

# the NWT LTER EDI id for observational data from the saddle

saddle_catch_sensntwk <- "210" # Saddle catchment sensor network data, 2017 - ongoing.
saddle_snow_depth_data <- "31" # Snow depth data for Saddle grid, 1992 - ongoing
saddle_productivity_data <- "16" # Aboveground net primary productivity data for Saddle (contains veg community classification for grid points) grid, 1992 - ongoing
saddle_sensntwk_veg <- "191" # Plot vegetation surveys at the Sensor Network, 2017 - ongoing
saddle_sensntwk_phenocam <- "192" # GCC data from Saddle sensor network phenocams, 2017 - ongoing
# Other possibly useful datasets:
# 211: Above-ground biomass for Sensor Node Array from 2017 to 2018, yearly
#      Has plots corresponding to sensor soil moisture node, also describes plot
#      biomass.
# 16:  Aboveground net primary productivity for saddle grid 
# not in edi yet: tvan soil temp/moisture from West tower. 
##############################################################################
# Helper functions - for downloading and loading data
##############################################################################

# Functions for downloading LTER Precip data are from Sarah Elmendorf's 
# utility_functions_all.R script
# https://github.com/NWTlter/long-term-trends/blob/master/utility_functions/utility_functions_all.R
# function to determine current version of data package on EDI
getCurrentVersion <- function(edi_id){
  require(magrittr)
  versions = readLines(paste0('https://pasta.lternet.edu/package/eml/knb-lter-nwt/', edi_id), 
                       warn = FALSE) %>%
    as.numeric() %>% (max)
  packageid = paste0('knb-lter-nwt.', edi_id, '.', versions)
  return(packageid)
}

#function to download the EML file from EDI
getEML <- function(packageid){
  require(magrittr)
  myurl <- paste0("https://portal.edirepository.org/nis/metadataviewer?packageid=",
                  packageid,
                  "&contentType=application/xml")
  #myeml<-xml2::download_html(myurl)%>%xml2::read_xml()%>%EML::read_eml()
  myeml <- xml2::read_xml(paste0("https://portal.edirepository.org/nis/metadataviewer?packageid=",
                                 packageid,
                                 "&contentType=application/xml")) %>% EML::read_eml()
}

# Function for downloading from EDI
download_EDI <- function(edi_id, dest_dir, getNewData = TRUE) {
  # This section heavily borrowed from Sarah Elmendorf's generic_timeseries_workflow.R script
  # https://github.com/NWTlter/long-term-trends/blob/master/plotting_scripts/generic_timeseries_workflow.R
  
  # Depends on getCurrentVersion() and getEML()
  
  packageid = getCurrentVersion(edi_id)
  
  # Check if files already exist; if they do, load the file paths and column classes
  # or warn the user if column classes can't be loaded
  if (any(grepl(packageid, list.files(dest_dir)) == TRUE)) {
    writeLines(paste0("Most recent package version ", 
                      packageid, " is already downloaded. Nothing to do."))
    # try to load the file path and column class information
    rds.fp <- list.files(dest_dir, pattern = ".RDS", full.names = TRUE)
    if (length(rds.fp) == 0) {
      warning(paste0("RDS file not found in download directory.",
                     " This means that no column class information is available,",
                     " so column classes will be set to NA."))
      csv_fp_colclass_list <- list(csv = list.files(dest_dir,
                            pattern = paste0(packageid, ".{1,}csv"), full.names = T),
                            colclasses = NA)
    } else {
      csv_fp_colclass_list <- readRDS(rds.fp)
    }
    
    return(csv_fp_colclass_list)
  # Check if a more recent version is available; if getNewData is set to false,
  # load column classes and file paths as above.
  } else if (getNewData == FALSE) {
    writeLines(paste0("A more recent version of the data (version ", 
                      packageid, ") is available. ",
                      "But since you have specified getNewData = FALSE, ",
                      "the latest version will not be downloaded."))
    if (length(rds.fp) == 0) {
      warning(paste0("RDS file not found in download directory.",
                     " This means that no column class information is available,",
                     " so column classes will be set to NA."))
      csv_fp_colclass_list <- list(csv = list.files(dest_dir, 
                                                    pattern = paste0(packageid, ".{1,}csv"), full.names = T),
                                   colclasses = NA)
    } else {
      csv_fp_colclass_list <- readRDS(rds.fp)
    }
    return(csv_fp_colclass_list)
  
  } 
  # Download the files if they need to be downloaded; extract the column classes from the 
  # attributes and save the file paths and column classes of the data in the dest_dir 
  else {
    
    writeLines(paste0("Downloading package ", packageid, " from EDI."))
    
    myeml = getEML(packageid)
    
    # Create output directory for data
    ifelse(!dir.exists(file.path(dest_dir)),
           dir.create(file.path(dest_dir)), FALSE)
    
    ### eml reading and downloading of csv
    if (is.null(names(myeml$dataset$dataTable))) {
      attributeList = lapply(myeml$dataset$dataTable, function(x){
        EML::get_attributes(x$attributeList)
      })
      names(attributeList) = lapply(myeml$dataset$dataTable, function(x){
        x$physical$objectName})
      if (getNewData) {
        #download all the datatables in the package
        csv_list <- list()
        csv_list <- lapply(myeml$dataset$dataTable, function(x){
          url_to_get = x$physical$distribution$online$url$url
          download.file(url_to_get,
                        destfile = paste0(dest_dir, "/",
                                          packageid, "_",
                                          x$physical$objectName),
                        method = "curl")
          output_csv <-  paste0(dest_dir, "/",
                                     packageid, "_",
                                     x$physical$objectName)
          return(output_csv)
          #print(output_csv_file)
        })
        output_csv_file <- csv_list
        
        # Get column classes
        colclasses <- list()
        colclasslist <- lapply(names(attributeList), 
                               function(x) {
                                 attributeList[[x]]$attributes = attributeList[[x]]$attributes %>%
                                   dplyr::mutate(
                                     col_classes = case_when(domain=='textDomain' ~'col_character()',
                                                           domain=='dateTimeDomain' ~'col_date()',
                                                           domain=='numericDomain' ~'col_number()',
                                                           domain=='enumeratedDomain' ~'character()')) %>%
                                   dplyr::mutate(
                                     colclasses=case_when(domain=='textDomain' ~'character',
                                                          domain=='dateTimeDomain' ~'Date',
                                                          domain=='numericDomain' ~'numeric',
                                                          domain=='enumeratedDomain' ~'character'))
                                 # replace "Date" with "Character" since R has a problem
                                 # with loading in Dates that aren't formatted exactly as
                                 # one wants
                                 attributeList[[x]]$attributes$colclasses <- gsub("Date",
                                                                                  "character",
                                                                                  attributeList[[x]]$attributes$colclasses)
                                 return(attributeList[[x]]$attributes$colclasses)
                               })
        names(colclasslist) <- names(attributeList)
        }
      
    }else{
      #if only one data table
      attributeList = list(EML::get_attributes(myeml$dataset$dataTable$attributeList))
      names(attributeList) = myeml$dataset$dataTable$physical$objectName
      if (getNewData) {
        url_to_get = myeml$dataset$dataTable$physical$distribution$online$url$url
        download.file(url_to_get,
                      destfile = paste0(dest_dir, "/",
                                        packageid, "_",
                                        myeml$dataset$dataTable$physical$objectName),
                      method = "curl")
        output_csv_file <- paste0(dest_dir, "/",
                                  packageid, "_",
                                  myeml$dataset$dataTable$physical$objectName)
      }
      # Get column classes for the data
      #map eml types to R col classes
      attributeList[[1]]$attributes <- attributeList[[1]]$attributes %>%
        dplyr::mutate(col_classes=case_when(domain=='textDomain' ~'col_character()',
                                            domain=='dateTimeDomain' ~'col_date()',
                                            domain=='numericDomain' ~'col_number()',
                                            domain=='enumeratedDomain' ~'character()')) %>%
        dplyr::mutate(colclasses=case_when(domain=='textDomain' ~'character',
                               domain=='dateTimeDomain' ~'Date',
                               domain=='numericDomain' ~'numeric',
                               domain=='enumeratedDomain' ~'character'))
      colclasslist <- attributeList[[1]]$attributes$colclasses
      
    }
    
    # Also save the full xml
    write_eml(myeml, file = paste0(dest_dir, "/", packageid, ".xml"))
    
    # Save the output list as R-object so column classes are preserved
    csv_fp_colclass_list <- list(csv = output_csv_file, colclasses = colclasslist)
    saveRDS(csv_fp_colclass_list, file = paste0(dest_dir, "/", packageid, 
                                                "_files_and_colclasses",
                                                ".RDS"))
    writeLines(paste0("Downloaded data can be found in: ", dest_dir))
    return(csv_fp_colclass_list)
  }
}

################################################################################
# Download Data
################################################################################
# Saddle sensor network
message(paste0("Downloading Saddle Catchment sensor network data, please cite: \n",
               "Morse, J. and Niwot Ridge LTER. 2020. Saddle catchment sensor network data, 2017- ongoing. ver 2. Environmental Data Initiative. https://doi.org/10.6073/pasta/9415ac5a669c11c6501612a94f90e04a (Accessed ",Sys.Date(), ")"))
saddle_catch_sensntwk_data_fp <- download_EDI(edi_id = saddle_catch_sensntwk,
                                              dest_dir = paste0(DirDnld,
                                                                "/saddle_sensorntwk_data"),
                                              getNewData = getNewData)
# Download sensor network community

# Download saddle grid snow_depth_data
message(paste0("Downloading Saddle Snow Depth data, please cite: \n",
               "Walker, S., J. Morse, and Niwot Ridge LTER. 2020. Snow depth data for Saddle grid, 1992 - ongoing ver 17. Environmental Data Initiative. https://doi.org/10.6073/pasta/8186d641539c37787495804b817e55ed (Accessed ",Sys.Date(), ")"))
saddle_snwdpt_data_fp <- download_EDI(edi_id = saddle_snow_depth_data, 
                                   dest_dir = paste0(DirDnld, "/saddle_snow_depth_data"),
                                   getNewData = getNewData)

# Download saddle grid productivity
message(paste0("Downloading Saddle Productivity data, please cite: \n",
               "Walker, M., J. Smith, H. Humphries, and Niwot Ridge LTER. 2019. Aboveground net primary productivity data for Saddle grid, 1992 - ongoing. ver 4. Environmental Data Initiative. https://doi.org/10.6073/pasta/34b6a7bbe47f9398ff7f5a748f90e838 (Accessed ",Sys.Date(), ")"))
saddle_prod_data_fp <- download_EDI(edi_id = saddle_productivity_data, 
                                   dest_dir = paste0(DirDnld, 
                                                     "/saddle_productivity_data"),
                                   getNewData = getNewData)

# Download saddle sensor network veg community
message(paste0("Downloading Saddle Productivity data, please cite: \n",
               "Elwood, K., W. Reed, and Niwot Ridge LTER. 2020. Plot vegetation surveys at the Sensor Network, 2017 to ongoing ver 2. Environmental Data Initiative. https://doi.org/10.6073/pasta/1b5e99d522f986c2244bf5a25e69d3f5 (Accessed ",Sys.Date(), ")"))
saddle_sensntwk_veg_data_fp <- download_EDI(edi_id = saddle_sensntwk_veg, 
                                    dest_dir = paste0(DirDnld, 
                                                      "/saddle_sensntwk_veg_data"),
                                    getNewData = getNewData)

# Download saddle sensor network phenocam data
message(paste0("Downloading Saddle Catchment phenocam data, please cite: \n",
               "Elwood, K., J. Smith, and Niwot Ridge LTER. 2021. Time-lapse camera (phenocam) imagery of Sensor Network plots from 2017 to ongoing. ver 2. Environmental Data Initiative. https://doi:10.6073/pasta/89e8189093392325ee139eccc6b2ff85 (Accessed ",Sys.Date(), ")"))
saddle_sensntwk_phenocam_data_fp <- download_EDI(edi_id = saddle_sensntwk_phenocam,
                                                 dest_dir = paste0(DirDnld,
                                                                   "/saddle_phenocam_data"),
                                                 getNewData = getNewData)

################################################################################
# Load Tvan flux data
################################################################################
# Both
tvan_comb <- read.table(file = tvan_data_fp, sep = "\t",
                        skip =  2, header = FALSE)
tvan_comb_names <- read.table(file = tvan_data_fp, sep = "\t",
                              header = TRUE, nrows = 1)
tvan_comb_units <- as.character(unname(unlist(tvan_comb_names[1,])))

colnames(tvan_comb) <- names(tvan_comb_names)
plot(tvan_comb$Tsoil)

# convert flux GPP (umol/m2/s to g/m2/s, as in CLM)
tvan_comb$GPP = tvan_comb$GPP * 1e-6 * 12.01
tvan_comb_units[1]  =  'gC m-2 s-1'
################################################################################
# Clean and format Tvan Flux data
################################################################################
tvan_comb_mod <- tvan_comb %>%
  mutate_all(list(~na_if(., -9999))) %>%
  mutate(timestamp = DateTime, 
         Hour = lubridate::hour(timestamp) + 
           lubridate::minute(timestamp)/60,
         date = lubridate::date(timestamp)) %>%
  mutate(DoY = yday(date), 
         Year = year(date),
         month = month(date)) %>%
  mutate(MonGroup = ifelse(month %in% c(12,1,2), "DJF",
                           ifelse(month %in% c(3,4,5), "MAM", 
                                  ifelse(month %in% c(6,7,8), "JJA", "SON")))) %>%
  #group_by(Year, DoY) %>%
  #mutate_at(all_of(c("NEE", "LE", "H", "Ustar", "Tair", "VPD", "rH", "VPD", "U",
  #                   "PRECTmms", "P", "FLDS", "Rg", "radNet", "Tsoil", "GPP")),
  #          list(daily_mean = mean), na.rm = TRUE) %>%
  select(date, timestamp, Year, DoY, Hour, everything())

# Get diurnal fluxes
# set the variables to use for diurnal fluxes
diurnal_flx_vars <- c("radNet", "H", "LE", "GPP")
tvan_comb_mod.diurnal_seasonal <- tvan_comb_mod %>% 
  select(-timestamp, -date) %>%
  select(Hour, DoY, Year, MonGroup, all_of(diurnal_flx_vars)) %>%
  group_by(MonGroup, Hour) %>%
  summarize_at(all_of(diurnal_flx_vars),
               list(houravg = mean, hoursd = sd), na.rm = TRUE) %>%
  mutate(ObsSim = "Obs") %>%
  mutate(veg_com = "FF")

# Get DoY fluxes
DoY_flx_vars <- c("GPP", "LE",'Tsoil')
tvan_comb_mod.daily <- tvan_comb_mod %>% 
  select(Hour, DoY, month, Year, all_of(DoY_flx_vars)) %>%
  # remove leap days and fix DoY
  filter(!(leap_year(Year) & DoY == 60)) %>%
  mutate(DoY = if_else(leap_year(Year) & (DoY > 59), 
                       DoY - 1, DoY)) %>% 
  group_by(DoY) %>%
  summarize_at(all_of(DoY_flx_vars),
               list(dailyavg = mean, dailysd = sd), na.rm = TRUE) %>%
  select(!starts_with("LE")) %>%
  mutate(ObsSim = "Obs") %>%
  mutate(veg_com = "FF")

plot(tvan_comb_mod.daily$Tsoil_dailyavg,type='l')

# Get July data
jul_30_min_tvan <- tvan_comb_mod %>% 
  select(-timestamp, -date) %>%
  select(Hour, DoY, Year, month, all_of(diurnal_flx_vars)) %>%
  filter(month == 7) %>%
  group_by(Hour) %>%
  summarize_at(all_of(diurnal_flx_vars),
               list(houravg = mean, hoursd = sd), na.rm = TRUE) %>%
  mutate(ObsSim = "Obs") %>%
  mutate(veg_com = "FF")

################################################################################
# Load Saddle Catchment Sensor Network Data
################################################################################
writeLines("Reading in saddle sensor network data...")

sad_sens_data_raw <- vector(length = length(saddle_catch_sensntwk_data_fp$csv),
                            mode = "list")
sad_sens_data_raw <- lapply(seq_along(saddle_catch_sensntwk_data_fp$csv), 
                            function(x) {
                              writeLines(paste0("Reading in ", 
                                                basename(saddle_catch_sensntwk_data_fp$csv[[x]])))
                              # replace date with chararacter
                              tmp_colclasses <- gsub("Date", "character", saddle_catch_sensntwk_data_fp$colclasses[[x]])
                              read.csv(saddle_catch_sensntwk_data_fp$csv[[x]],
                                       header = T, sep = ",", quot = '"', 
                                       as.is = TRUE, na.strings = c("NA", "NaN", ""),
                                       colClasses = tmp_colclasses)
                            })
names(sad_sens_data_raw) <- basename(unlist(saddle_catch_sensntwk_data_fp$csv))

sad_sens_data_all <- bind_rows(sad_sens_data_raw)

#writeLines("Reading in saddle sensor veg data...")


################################################################################
# Gather vegetation communities for Saddle Network plots
################################################################################
# Using Kelsey Elwood Carter's Community characterizations from her master's thesis
sensor_plot_com <- data.frame(plot = c(9,10,14, 16, 20, 21, 11, 15, 6,7,8,12,13,17,19),
                              veg_com_long = c(rep("Dry Meadow 1", 3), 
                                               rep("Dry Meadow 2", 3),
                                               rep("Dry Meadow 3", 2),
                                               rep("Moist Meadow", 4),
                                               rep("Wet Meadow", 2),
                                               "Subalpine")) %>%
  mutate(veg_com = ifelse(grepl("Dry", veg_com_long), "DM",
                   ifelse(grepl("Moist", veg_com_long), "MM",
                   ifelse(grepl("Wet", veg_com_long), "WM", 
                   ifelse(grepl("Subalpine", veg_com_long), "SA", NA)))),
         plot = as.character(plot))

  
################################################################################
# Handle Saddle Network Soil Moisture and Temperature data
################################################################################
# Filter out questionable data from sensor network and collapse 10-minute readings
# into 30-minute readings. Also categorize the vegetation communities for each
# site
sad_sens_10min <- sad_sens_data_all %>% 
  mutate(timestamp = as.POSIXct(date, format = "%Y-%m-%d %H:%M:%OS", tz = "MST"),
         date = as.Date(timestamp)) %>% 
  select(sensornode, timestamp, date, contains("soil")) %>% 
  #filter(sensornode == 6) %>%
  # set flagged values to NA so they won't be used in averages
  mutate(soiltemp_5cm_avg = ifelse(!is.na(flag_soiltemp_5cm_avg), NA, soiltemp_5cm_avg),
         soiltemp_30cm_avg = ifelse(!is.na(flag_soiltemp_30cm_avg), NA, soiltemp_30cm_avg),
         soilmoisture_a_5cm_avg = ifelse(!is.na(flag_soilmoisture_a_5cm_avg), NA,
                                         soilmoisture_a_5cm_avg),
         soilmoisture_a_30cm_avg = ifelse(!is.na(flag_soilmoisture_a_30cm_avg), NA,
                                         soilmoisture_a_30cm_avg),
         soilmoisture_b_5cm_avg = ifelse(!is.na(flag_soilmoisture_b_5cm_avg), NA,
                                         soilmoisture_b_5cm_avg),
         soilmoisture_b_30cm_avg = ifelse(!is.na(flag_soilmoisture_b_30cm_avg), NA,
                                          soilmoisture_b_30cm_avg),
         soilmoisture_c_5cm_avg = ifelse(!is.na(flag_soilmoisture_c_5cm_avg), NA,
                                         soilmoisture_c_5cm_avg),
         soilmoisture_c_30cm_avg = ifelse(!is.na(flag_soilmoisture_c_30cm_avg), NA,
                                          soilmoisture_c_30cm_avg)) %>% 
  # Remove flag columns
  select(-contains("flag")) %>% 
  # Remove soil moisture data when temperature <0 (frozen water messes with sensors)
  mutate(soilmoisture_a_5cm_avg = ifelse(soiltemp_5cm_avg <= 0, NA,
                                         soilmoisture_a_5cm_avg),
         soilmoisture_a_30cm_avg = ifelse(soiltemp_30cm_avg <= 0, NA,
                                         soilmoisture_a_30cm_avg),
         soilmoisture_b_5cm_avg = ifelse(soiltemp_5cm_avg <= 0, NA,
                                         soilmoisture_b_5cm_avg),
         soilmoisture_b_30cm_avg = ifelse(soiltemp_30cm_avg <= 0, NA,
                                         soilmoisture_b_30cm_avg),
         soilmoisture_c_5cm_avg = ifelse(soiltemp_5cm_avg <= 0, NA,
                                         soilmoisture_c_5cm_avg),
         soilmoisture_c_30cm_avg = ifelse(soiltemp_30cm_avg <= 0, NA,
                                         soilmoisture_c_30cm_avg)) %>%
  # Group times by half-hour so half-hourly averages can be taken
  mutate(Time = gsub(".{4}-.{2}-.{2} ", "", timestamp),
           cleanTime =
           strsplit(Time, ":") %>%
           sapply(function(x){
             x <- as.numeric(x)
             x[1] + x[2]/60 + x[3]/(60*60)
           }),
         decimalTime = floor(cleanTime * 2)/2)

writeLines(paste0("Collapsing 10-minute soil sensor data into 30-minute chunks, \n", 
           "this may take a while..."))

# NOTE: this could probably be made more efficient if handled one file at a time. And then
# joining the 30-minute data together after each is combined
sad_sens_soilmoist_temp <- sad_sens_10min %>%
  # Get half-hourly averages
  group_by(date, decimalTime, sensornode) %>%
  mutate(across(contains("soil"), list(~mean(., na.rm = TRUE)), 
                .names = "mean_{col}")) %>%
  ungroup() %>%
  select(sensornode, date, decimalTime, contains("mean_")) %>%
  unique() %>%
  # Join with vegetation classifications
  left_join(sensor_plot_com, by = c("sensornode" = "plot")) %>%
  # Remove sub alpine, and non-characterized vegetation communities
  filter(!is.na(veg_com)) %>%
  filter(veg_com != "SA")


# Subset the data so that the "a" sensor is used, unless that sensor is NA, then
# preferentially use "b", and then "c"
sad_sensnet_soil <- sad_sens_soilmoist_temp %>% 
  # if a is NA choose b
  mutate(soilmoisture_5cm_sensor_letter = ifelse(is.na(mean_soilmoisture_a_5cm_avg),
                                                     "b", "a"),
         soilmoisture_5cm_avg = ifelse(is.na(mean_soilmoisture_a_5cm_avg),
                                       mean_soilmoisture_b_5cm_avg, 
                                       mean_soilmoisture_a_5cm_avg),
         soilmoisture_30cm_sensor_letter = ifelse(is.na(mean_soilmoisture_a_30cm_avg),
                                                 "b", "a"),
         soilmoisture_30cm_avg = ifelse(is.na(mean_soilmoisture_a_30cm_avg),
                                       mean_soilmoisture_b_30cm_avg, 
                                       mean_soilmoisture_a_30cm_avg)) %>%
  # if b is also NA choose c
  mutate(soilmoisture_5cm_sensor_letter = ifelse(is.na(soilmoisture_5cm_avg),
                                                 "c", soilmoisture_5cm_sensor_letter),
         soilmoisture_5cm_avg = ifelse(is.na(soilmoisture_5cm_avg),
                                       mean_soilmoisture_c_5cm_avg, 
                                       soilmoisture_5cm_avg),
         soilmoisture_30cm_sensor_letter = ifelse(is.na(soilmoisture_30cm_avg),
                                                 "c", soilmoisture_30cm_sensor_letter),
         soilmoisture_30cm_avg = ifelse(is.na(soilmoisture_30cm_avg),
                                       mean_soilmoisture_c_30cm_avg, 
                                       soilmoisture_30cm_avg)) %>%
  # if c is also NA, make everything NA
  mutate(soilmoisture_5cm_sensor_letter = ifelse(is.na(soilmoisture_5cm_avg),
                                                 NA, soilmoisture_5cm_sensor_letter),
         soilmoisture_5cm_avg = ifelse(is.na(soilmoisture_5cm_avg),
                                       NA, 
                                       soilmoisture_5cm_avg),
         soilmoisture_30cm_sensor_letter = ifelse(is.na(soilmoisture_30cm_avg),
                                                 NA, soilmoisture_30cm_sensor_letter),
         soilmoisture_30cm_avg = ifelse(is.na(soilmoisture_30cm_avg),
                                       NA, 
                                       soilmoisture_30cm_avg)) %>% 
  select(sensornode, date, decimalTime, mean_soiltemp_5cm_avg, mean_soiltemp_30cm_avg,
         soilmoisture_5cm_avg, soilmoisture_30cm_avg, veg_com, 
         soilmoisture_5cm_sensor_letter,
         soilmoisture_30cm_sensor_letter) %>%
  # rename to be generic enough to match Tvan soil columns
  dplyr::rename(soiltemp_upper_avg = mean_soiltemp_5cm_avg, 
         soiltemp_lower_avg = mean_soiltemp_30cm_avg,
         soilmoisture_upper_avg = soilmoisture_5cm_avg, 
         soilmoisture_lower_avg = soilmoisture_30cm_avg,
         sensnet_soilmoisture_5cm_sensor_letter = soilmoisture_5cm_sensor_letter,
         sensnet_soilmoisture_30cm_sensor_letter = soilmoisture_30cm_sensor_letter,
         Hour = decimalTime,
         plot = sensornode) %>%
  mutate(data_set = "Saddle_sensor_network_EDI_210_5cm_30cm_moisttemp_probes",
         plot = as.character(plot),
         upper_sensor_depth_cm = 5,
         lower_sensor_depth_cm = 30) %>%
  mutate(soilmoisture_upper_avg = soilmoisture_upper_avg * 100,
         soilmoisture_lower_avg = soilmoisture_lower_avg * 100)


# Plotting soil moisture
# sad_sens_soilmoist_temp.nested <- sad_sens_soilmoist_temp %>%
#   mutate(timestamp = as.Date(date, origin = paste0(date, " 00:00:00")) +
#            lubridate::minutes(decimalTime * 60)) %>%
#   #select(sensornode, date, decimalTime, timestamp) %>%
#   pivot_longer(contains("soil"), names_to = "Variable", values_to = "Value") %>%
#   group_by(veg_com) %>%
#   nest()

# sad_sens_soilmoist_temp.plot <- sad_sens_soilmoist_temp.nested %>%
#   mutate(plot = map2(data, veg_com, ~ ggplot(data = .x,
#                                              aes(x = timestamp, y = Value)) +
#                        ggtitle(glue::glue("Vegetation Community: {.y}")) +
#                        geom_point(alpha = 0.3) +
#                        facet_wrap(~Variable, scales = "free_y", ncol = 2)
#                      )
#          )
# 
# pdf("~/Downloads/test.pdf")
# sad_sens_soilmoist_temp.plot$plot
# dev.off()


################################################################################
# Tvan soil moisture data
################################################################################
# Read in Tvan soil moisture and temperature data
tvan_soil <- read.table(file = tvan_data_soil, sep = ",",
                        skip =  2, header = FALSE)
tvan_soil_names <- read.table(file = tvan_data_soil, sep = ",",
                              header = TRUE, nrows = 1)
tvan_soil_units <- as.character(unname(unlist(tvan_soil_names[1,])))

colnames(tvan_soil) <- names(tvan_soil_names)

# Fix time, rename to match generic names of sensor network soil data, add informational
# columns about the data's origin, and vegetation community
tvan_soil_mod <- tvan_soil %>%
  select(time, wc10, wc30, soil_temp, tc30, G) %>%
  mutate(timestamp = time,
         Hour = lubridate::hour(timestamp) +
           lubridate::minute(timestamp)/60,
         date = lubridate::date(timestamp)) %>%
  select(-time, -timestamp) %>%
  rename(soiltemp_upper_avg = soil_temp,
         soiltemp_lower_avg = tc30,
         soilmoisture_upper_avg = wc10,
         soilmoisture_lower_avg = wc30) %>%
  mutate(upper_sensor_depth_cm = 10,
         lower_sensor_depth_cm = 30) %>%
  mutate(veg_com = "FF",
         plot = "Tvan_West",
         data_set = "Tvan_West_Tower_10cm_30cm_moisttemp_probes")


plot(tvan_soil_mod$date, tvan_soil_mod$soiltemp_upper_avg,pch='.')
#ggplot(tvan_soil_mod, aes(x = DoY, y = soiltemp_upper_avg)) + geom_point() 


################################################################################
# Combine Tvan and Sensor Network Soil and Moisture data
################################################################################
soilmoist_temp_comb <- full_join(sad_sensnet_soil, tvan_soil_mod, 
                                 by = c("date", "Hour", 
                                        "veg_com", "plot", "data_set", 
                                        "soiltemp_upper_avg",
                                        "soiltemp_lower_avg", 
                                        "soilmoisture_upper_avg",
                                        "soilmoisture_lower_avg",
                                        "upper_sensor_depth_cm",
                                        "lower_sensor_depth_cm"))

# Summarize the data by hour
soilmoist_temp_comb_hrly <- soilmoist_temp_comb %>% 
  group_by(Hour, veg_com) %>%
  mutate_at(all_of(c("soiltemp_upper_avg",
                      "soiltemp_lower_avg",
                      "soilmoisture_upper_avg",
                      "soilmoisture_lower_avg")),
            list(~mean(., na.rm = TRUE), ~sd(., na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(soilmoisture_upper_avg = ifelse(soiltemp_upper_avg < 0, NA, 
                                         soilmoisture_upper_avg),
         soilmoisture_lower_avg = ifelse(soiltemp_lower_avg < 0, NA, 
                                         soilmoisture_lower_avg)) %>%
  mutate(data_information = paste0("data = ", data_set, " | ", 
                                   "sensnet_5cm_letter = ", 
                                   sensnet_soilmoisture_5cm_sensor_letter, " | ", 
                                   "sensnet_30cm_letter = ", 
                                   sensnet_soilmoisture_30cm_sensor_letter, " | ",
                                   "upr_sens_depth = ", upper_sensor_depth_cm, " | ",
                                   "lwr_sens_depth = ", lower_sensor_depth_cm)) %>%
  select(Hour, ends_with("_avg_mean"), ends_with("_avg_sd"), veg_com, 
         data_information) %>%
  unique()
  
# Summarize the data by daily averages
soilmoist_temp_comb_daily <- soilmoist_temp_comb %>% 
  mutate(DoY = yday(date),
         month = month(date),
         year = year(date)) %>%
  # remove leap days and fix DoY
  filter(!(leap_year(year) & DoY == 60)) %>%
  mutate(DoY = if_else(leap_year(year) & (DoY > 59), 
                       DoY - 1, DoY)) %>%
  group_by(DoY, veg_com) %>%
  mutate(across(all_of(c("soiltemp_upper_avg",
                     "soiltemp_lower_avg",
                     "soilmoisture_upper_avg",
                     "soilmoisture_lower_avg")),
                .fns = list(dailyavg = ~mean(., na.rm = TRUE), 
                     dailysd = ~sd(., na.rm = TRUE)))) %>%
  ungroup() %>%
  mutate(soilmoisture_upper_avg = ifelse(soiltemp_upper_avg < 0, NA, 
                                         soilmoisture_upper_avg),
         soilmoisture_lower_avg = ifelse(soiltemp_lower_avg < 0, NA, 
                                         soilmoisture_lower_avg)) %>%
  # construct a data information column with information about which dataset and
  # sensor the data came from
  mutate(data_information = paste0("data = ", data_set, " | ", 
                                   "sensnet_5cm_letter = ", 
                                   sensnet_soilmoisture_5cm_sensor_letter, " | ", 
                                   "sensnet_30cm_letter = ", 
                                   sensnet_soilmoisture_30cm_sensor_letter, " | ",
                                   "upr_sens_depth = ", upper_sensor_depth_cm, " | ",
                                   "lwr_sens_depth = ", lower_sensor_depth_cm)) %>%
  select(DoY, month, ends_with("_avg_dailyavg"), ends_with("_avg_dailysd"), 
         veg_com, data_information) %>%
  unique()

names(soilmoist_temp_comb_daily)
ggplot(soilmoist_temp_comb_daily, aes(x = DoY)) +
  geom_line(aes(y = soiltemp_upper_avg_dailyavg, color = veg_com))


################################################################################
# Load in data from Saddle Grid (snow depth and productivity)
################################################################################
writeLines("Reading in saddle grid snow depth data...")
# Daily data
sad_snw <- read.csv(saddle_snwdpt_data_fp$csv,
                    header = T, sep = ",", quot = '"')

writeLines("Reading in saddle productivity data...")

sad_prod <- read.csv(saddle_prod_data_fp$csv[[1]],
                     header = T, sep = ",", quot = '"')

################################################################################
# Handle Saddle Grid Snow-depth data
################################################################################
# Get saddle grid point vegetation community characterizations
sad_grid_veg_com <- sad_prod %>%
  select(grid_pt, veg_class) %>%
  dplyr::rename(veg_com = veg_class) %>%
  unique()

# Merge saddle snow depth measurements with the saddle vegetation characterizations
# Snow measured in cm at stake and m in the model (SNOW_DEPTH)
sad_snw_mod <- sad_snw %>%
  left_join(sad_grid_veg_com, by = c("point_ID" = "grid_pt")) %>%
  #select(point_ID, veg_class) %>%
  mutate(veg_class = ifelse(is.na(veg_com), "not available", veg_com)) %>%
  filter(!(veg_com == "not available")) %>%
  filter(veg_com %in% c("DM", "FF", "MM", "SB", "WM")) %>%
  mutate(snow_depth = as.numeric(mean_depth),
         date = as.Date(date, format = "%Y-%m-%d"),
         DoY = lubridate::yday(date),
         Year = lubridate::year(date)) %>%
  group_by(veg_com, DoY) %>%
  mutate(snow_depth_dailyavg = mean(snow_depth, na.rm=TRUE),
         snow_depth_dailysd = sd(snow_depth, na.rm = TRUE)) %>%
  select(date, DoY, Year, point_ID, snow_depth, snow_depth_dailyavg,
         snow_depth_dailysd, veg_com) %>% 
  ungroup() %>%
  mutate(data_information = "Saddle_grid_snow_depth_EDI_31")

# Get DoY averages
sad_snw_daily <- sad_snw_mod %>%
  select(DoY, snow_depth_dailyavg, snow_depth_dailysd,
         veg_com, data_information) %>%
  unique() %>%
  dplyr::rename(snow_depth_data_information = data_information)

# Average the snow depth across plots of the same vegetation community, at 
# each date
sad_snw_forc_yrs <- sad_snw_mod %>%
  filter(Year >= 2007) %>%
  group_by(date, veg_com) %>%
  mutate(avg_date_depth = mean(snow_depth, na.rm = TRUE),
         sd_date_depth = sd(snow_depth, na.rm = TRUE)) %>%
  ungroup() %>%
  select(date, DoY, Year, avg_date_depth, sd_date_depth, veg_com,
         data_information) %>%
  unique()



# plot the measurements and doy averages for each community
# sad_snw_mod %>%
#   filter(snow_depth > 0) %>%
#   ggplot(aes(x = as.Date("2000-01-01", format = "%Y-%m-%d") + 
#                (DoY - 1))) +
#   #geom_point(aes(y = snow_depth), alpha = 0.03) +
#   geom_line(data = sad_snw_mod %>%
#               select(DoY, doy_avg_depth, veg_com) %>%
#               unique(),
#             aes(y = doy_avg_depth), color = "red") +
#   geom_ribbon(aes(ymin = doy_avg_depth - doy_sd_depth, 
#                   ymax = doy_avg_depth + doy_sd_depth),
#               alpha = 0.3) +
#   facet_wrap(~veg_com) +
#   scale_x_date(date_labels = "%b", date_breaks = "1 month")
# 
# sad_snw_mod %>%
#   filter(Year >= 2007) %>%
#   group_by(date, veg_com) %>%
#   mutate(avg_date_depth = mean(snow_depth, na.rm = TRUE),
#          sd_date_depth = sd(snow_depth, na.rm = TRUE)) %>%
#   ungroup() %>%
#   ggplot(aes(x = date)) +
#   geom_ribbon(aes(ymin = avg_date_depth - sd_date_depth,
#                   ymax = avg_date_depth + sd_date_depth),
#               alpha = 0.5) +
#   geom_line(aes(y = avg_date_depth)) +
#   facet_wrap(~veg_com, ncol = 1) +
#   scale_x_date(date_labels = "%b-%Y", date_breaks = "1 year") +
#   theme(axis.text.x = element_text(angle = 45, hjust =  1))

  
################################################################################
# Handle Saddle Grid Productivity data
################################################################################
# Yearly Saddle grid Productivity data (gC/m^2)
# CLM needs it in gC/m^2/s
sad_prod_mod <- sad_prod %>% 
  dplyr::rename(veg_com = veg_class) %>%
  select(year, grid_pt, veg_com, NPP, subsample) %>%
  #mutate(row = row_number()) %>%
  # Separate by subsamples
  pivot_wider(names_from = matches("subsample"),
              names_prefix = "subsample_",
              values_from = matches("NPP")) %>%
  # average subsamples
  mutate(NPP = rowMeans(select(., starts_with("subsample_")),
                    na.rm = TRUE))
sad_prod_mod_ann <- sad_prod_mod %>%
  group_by(year, veg_com) %>%
  mutate(mean_NPP = mean(NPP, na.rm = TRUE),
         sd_NPP = sd(NPP, na.rm = TRUE))
  


# sad_prod_mod %>%
#   ggplot(aes(x = veg_com, y = NPP)) +
#   geom_boxplot(fill = NA) +
#   geom_point(position = position_jitter(width = rel(0.3)))

# Get relevant meadow measurements: -- use purr to put in list; filter out tundra shrub and
# snow fence first
# sad_prod_FF <- sad_prod_sub %>%
#   filter(veg_class == "FF")
# sad_prod_SB <- sad_prod_sub %>%
#   filter(veg_class == "SB")
# sad_prod_MM <- sad_prod_sub %>%
#   filter(veg_class == "MM")
# sad_prod_WM <- sad_prod_sub %>%
#   filter(veg_class == "WM")
# sad_prod_DM <- sad_prod_sub %>%
#   filter(veg_class == "DM")


################################################################################
# Load Saddle Catchment Phenocam Data
################################################################################
writeLines("Reading in sensor network phenocam data...")

phenocam_data_raw <- vector(length = length(saddle_sensntwk_phenocam_data_fp$csv),
                                            mode = "list")
phenocam_data_raw <- lapply(seq_along(saddle_sensntwk_phenocam_data_fp$csv), 
                                            function(x) {
                                              writeLines(paste0("Reading in ", 
                                                                basename(saddle_sensntwk_phenocam_data_fp$csv[[x]])))
                                              # replace date with chararacter
                                              tmp_colclasses <- gsub("Date", "character", saddle_sensntwk_phenocam_data_fp$colclasses[[x]])
                                              read.csv(saddle_sensntwk_phenocam_data_fp$csv[[x]],
                                                       header = T, sep = ",", quot = '"', 
                                                       as.is = TRUE, na.strings = c("NA", "NaN", ""),
                                                       colClasses = tmp_colclasses)
                                            })
names(phenocam_data_raw) <- basename(unlist(saddle_sensntwk_phenocam_data_fp$csv))

phenocam_data_filtered <- phenocam_data_raw[[2]]

# Joining veg classifications with phenocam data
phenocam_data_filtered <- left_join(phenocam_data_filtered, sensor_plot_com, by = c("node" = "plot")) %>%
  # Remove sub alpine, and non-characterized vegetation communities
  filter(!is.na(veg_com)) %>%
  filter(veg_com != "SA")

################################################################################
# Reformat data
################################################################################
# Reformatting several data frames to better match with simulation data

# Data frame 1: 
# Rename flux variables to match tvan
# Half-hourly fluxes from Tvan; Comparable to the fell-field
# Variables: FSH (tvan), RN (tvan), LE (tvan), GPP (tvan)
tvan_comb_mod.diurnal_seasonal <- tvan_comb_mod.diurnal_seasonal %>%
  dplyr::rename(RNET_houravg = radNet_houravg,
         RNET_hoursd = radNet_hoursd,
         FSH_houravg = H_houravg,
         FSH_hoursd = H_hoursd,
         EFLX_LH_TOT_houravg = LE_houravg,
         EFLX_LH_TOT_hoursd = LE_hoursd)


# July flux summary
jul_30_min_tvan <- jul_30_min_tvan %>%
  dplyr::rename(RNET_houravg = radNet_houravg,
         RNET_hoursd = radNet_hoursd,
         FSH_houravg = H_houravg,
         FSH_hoursd = H_hoursd,
         EFLX_LH_TOT_houravg = LE_houravg,
         EFLX_LH_TOT_hoursd = LE_hoursd)

# Data frame 2: 
# Daily averages for each vegetation community
# Variables: GPP (tvan), SoilTemp (tvan/sensor network),
# Soil Moisture (tvan/sensor network), snow depth (saddle grid), 
daily_soilmoisttemp_gpp_snwdp <- soilmoist_temp_comb_daily %>% 
  dplyr::rename(soilmoisture_data_info = data_information) %>%
  # join with tvan GPP data
  left_join(tvan_comb_mod.daily, by = c("DoY", "veg_com")) %>%
  # join with snow depth data
  left_join(sad_snw_daily, by = c("DoY", "veg_com")) %>%
  mutate(ObsSim = "Obs")
  


################################################################################
# Combine and write out
################################################################################
# For each time-series of data, write out units and data
# Write out halfhourly fluxes:
writeLines("Writing out diurnal, daily, and annual data.")

# Diurnal-seasonal data
write.table(tvan_comb_mod.diurnal_seasonal, 
            file = paste0(DirOutBase, "/Diurnal_seasonal_summaries_", "tvan_flux.txt"),
            row.names = FALSE, sep = "\t")

# Diurnal-seasonal data
write.table(jul_30_min_tvan, 
            file = paste0(DirOutBase, "/July_flux_summary_", "tvan_flux.txt"),
            row.names = FALSE, sep = "\t")

# DoY data 
write.table(daily_soilmoisttemp_gpp_snwdp, 
            file = paste0(DirOutBase, 
                          "/Daily_soilmoisture_soiltemp_gpp_snwdpth_summaries.txt"),
            row.names = FALSE, sep = "\t")


# Annual data
write.table(sad_prod_mod_ann, 
            file = paste0(DirOutBase, 
                          "/annual_saddle_grid_NPP_summaries.txt"),
            row.names = FALSE, sep = "\t")

# Unsummarized data
writeLines("Writing out data that has not been summarized by time.")

# Saddle sensor network soil data
write.table(sad_sensnet_soil, 
            file = paste0(DirOutBase, 
                          "/sensor_network_soil_data_30_min.txt"),
            row.names = FALSE, sep = "\t")

# Tvan soil data
write.table(tvan_soil_mod,
            file = paste0(DirOutBase, 
                          "/tvan_soil_data_30_min.txt"),
            row.names = FALSE, sep = "\t")

# Snow depth data
write.table(sad_snw_forc_yrs,
            file = paste0(DirOutBase, 
                          "/saddle_grid_snow_depth_data_biweekly.txt"),
            row.names = FALSE, sep = "\t")

# Productivity
write.table(sad_prod_mod,
            file = paste0(DirOutBase, 
                          "/saddle_grid_productivity_data.txt"),
            row.names = FALSE, sep = "\t")

# Phenocam data
write.table(phenocam_data_filtered,
            file = paste0(DirOutBase, 
                          "/sensor_network_phenocam_data.txt"),
            row.names = FALSE, sep = "\t")

print('--- finished  with script ---')
