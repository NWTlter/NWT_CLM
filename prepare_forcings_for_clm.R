##############################################################################################
#' title Workflow to NCAR CLM data set

#' author
#' Hannah Holland-Moritz (hhollandmoritz AT gmail.com), based on script by David Durden (eddy4R.info AT gmail.com)
#' 
#' description 
#' Workflow for collating NIWOT LTER data, gap-filling, and packaging in NCAR CLM netcdf format.

# Modified from David Durden's flow.api.clm.R script for NEON data
# changelog and author contributions / copyrights
# David Durden (2019-07-05)
#   original creation
# David Durden (2020-05-31)
# Updating to use neonUtilities for all data retrieval from API
##############################################################################

##############################################################################
# Dependencies
##############################################################################

#Call the R HDF5 Library
packReq <- c("rhdf5","REddyProc", "ncdf4","devtools","magrittr","EML", "dplyr",
             "ggplot2", "purrr", "tidyr", "lubridate","RCurl", "httr", "jsonlite")

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
#### Ploting options ####
# Should plots be made of gap-filled data?
makeplots <- TRUE # FALSE

#### Output Options ####
# Base directory for all files
DirBase <- "~/Desktop/Working_files/Niwot/"
# Base directory for output
DirOutBase <- paste0(DirBase,"CLM/data")

#### Download and input options ####
# Directory to download precipitation and radidation data to
DirDnld = paste0(DirBase,"lter_flux")

# Should a newer version of precip data be automatically 
# downloaded if one is available?
getNewData = TRUE
# Ameriflux username
# NOTE: you cannot download Ameriflux data without a valid username
# to create an account, visit the Ameriflux website: https://ameriflux.lbl.gov/
# Please also read their data-use policy, by downloading their data you are agreeing 
# to follow it. The policy can be found here: https://ameriflux.lbl.gov/data/data-policy/
amf_usr <- "wwieder" # CHANGE ME

#### Tower Use Options ####
# What tvan tower should be used?
tower <- "Both" # Options are "East", "West", or "Both"
# if "Both" the one tower will be used to gapfill the other tower
# basetower provides which tower is the baseline that will be filled
# with the other tower. Currently the East tower record is more complete
# and has fewer gaps and errors, so it is being used as the basetower.
basetower <- "East" # West


#### Tvan data location ####
# Only necessary to set the location of the tower that you are processing, or 
# both, if tower = "Both"
# The data should be formatted with ReddyProc file format.
# Briefly the file should be formated as follows: the file should be 
# tab-delimited with the first row specifying the name of the variable 
# and the second specifying the units of that variable. The columns should have names 
# and units that follow the guidelines below:
# Column formating guidelines for Tvan data
# (optional indicates a column is not necessary for producing the final netcdf, 
# it includes variables that are necessary for CLM, and also variables that are
# necessary for ReddyProc gapfilling of the data in preparation for CLM).
# | Column Name |        Column Description        | Units          | Optional? |
# | ----------- | -------------------------------- | -------------- | --------- |
# | NEE         | Net ecosystem exchange           | umol m^-2 s^-1 | Yes       |
# | LE          | Latent heat flux                 | W m^-2         | No        |
# | H           | Sensible heat flux               | W m^-2         | No        |
# | Ustar       | Friction velocity                | m s^-1         | Yes       |
# | Tair        | Air temperature                  | degC           | No        |
# | VPD         | Vapor pressure density           | kPa            | No        |
# | rH          | relative humidity                | %              | No        |
# | U           | Wind speed                       | m s^-1         | No        |
# | P           | Atmospheric pressure             | kPa            | No        |
# | Tsoil       | Soil temperature                 | degC           | Yes       |
# | Year        | Year                             | -              | No        |
# | DoY         | The day of year (1-365/366)      | -              | No        |
# | Hour        | Decimal hour of the day (0.5-24) | -              | No        |
# The location of the east tvan data filepath, use "", if tower = "West"
DirIN = paste0(DirBase,"Tvan_out_new/supp_filtering/")
east_data_fp <- paste0(DirIN,"tvan_East_2007-05-10_00-30-00_to_2021-03-02_flux_P_reddyproc_cleaned.txt")
# The location of the west tvan data filepath, use "", if tower = "East"
west_data_fp <- paste0(DirIN,"tvan_West_2007-05-10_00-30-00_to_2021-03-02_flux_P_reddyproc_cleaned.txt")

#### Simulated Runoff Option ####
# WARNING THIS FEATURE IS UNTESTED; CHANGE AT YOUR OWN RISK
# The user can provide a data file from a simulated Moist Meadow run that
# contains two columns, a timestamp column (every timestamp represents the
# state at the *end* of the 30 minute sampling period) called "time",
# and a column containing the QRUNOFF amounts in mm/s from a Moist Meadow 

# simulation. If provided, this data will be added to the Wet meadow 
# precipitation. If not provided, wet meadow precipitation will be 75% of 
# observed precipitation.
# As done in Wieder et al. 2017, JGR-B. doi:10.1002/2016JG003704.

# Provide a character string specifying the location of the simulated runoff data
# if NA, no simulated runoff will be used
simulated_runoff_fp <- paste0(DirIN,'QRUNOFF_clm50bgc_NWT_mm_newPHS_lowSLA.csv')


##############################################################################
# Static workflow parameters - these are unlikely to change
##############################################################################

#Append the site to the base output directory
DirOut <- paste0(DirOutBase, "/", "data")
plots_dir <- paste0(DirOutBase, "/plots")

# Check if directory exists and create if not
if (!dir.exists(DirOut)) dir.create(DirOut, recursive = TRUE)
if (!dir.exists(DirDnld)) dir.create(DirDnld, recursive = TRUE)
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)


# the EDI id for precip data from the saddle and C1 weather stations
saddle_precip_data <- "416" # NWT LTER EDI id

# Lat/long coords - shouldn't need to change unless modified in surface
# dataset lat/long
latSite <- 40.05 # should match the lat of the surface dataset
lonSite <- 360 - 254.42 #  should match the long of the surface dataset

# Should simulated runoff mode be activated?
if (is.na(simulated_runoff_fp)) {
  simulated_runoff_present <- FALSE
  writeLines(paste0("No simulated runoff file supplied. Wet meadow precipitation",
                    "  will be calculated without any added runoff."))
} else {
  simulated_runoff_present <- TRUE
  writeLines(paste0("You have supplied the following simulated runoff file: \n",
                    simulated_runoff_fp,
                    "\nIt will be added when wet meadow precipitation", 
                    " is calculated."))
}


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
  myurl<-paste0("https://portal.edirepository.org/nis/metadataviewer?packageid=",
                packageid,
                "&contentType=application/xml")
  #myeml<-xml2::download_html(myurl)%>%xml2::read_xml()%>%EML::read_eml()
  myeml<-xml2::read_xml(paste0("https://portal.edirepository.org/nis/metadataviewer?packageid=",

                               packageid,
                               "&contentType=application/xml")) %>% EML::read_eml()
}

# Function for downloading from EDI
download_EDI <- function(edi_id, dest_dir, getNewData = TRUE) {
  # This section heavily borrowed from Sarah Elmendorf's generic_timeseries_workflow.R script
  # https://github.com/NWTlter/long-term-trends/blob/master/plotting_scripts/generic_timeseries_workflow.R
  
  # Depends on getCurrentVersion() and getEML()
  
  packageid = getCurrentVersion(edi_id)
  
  if (any(grepl(packageid, list.files(dest_dir)) == TRUE)) {
    writeLines(paste0("Most recent package version ", 
                      packageid, " is already downloaded. Nothing to do."))
    return(list.files(dest_dir, pattern = paste0(packageid, ".{1,}csv"), full.names = T))
  } else if (getNewData == FALSE) {
    writeLines(paste0("A more recent version of the data (version ", 
                      packageid, ") is available. ",
                      "But since you have specified getNewData = FALSE, ",
                      "the latest version will not be downloaded."))
    return(list.files(dest_dir, pattern = paste0(".{1,}csv"), full.names = T))
  } else {
    
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
                                        myeml$dataset$dataTable$physical$objectName),
                        method = "curl")
          output_csv_file <-  paste0(dest_dir, "/",
                                     packageid, "_",
                                     myeml$dataset$dataTable$physical$objectName)
        })
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
    }
    
    # Also save the full xml
    write_eml(myeml, file = paste0(dest_dir, "/", packageid, ".xml"))
    
    writeLines(paste0("Downloaded data can be found in: ", dest_dir))
    return(output_csv_file)
  }
}


# Function for downloading USCRN precip
download_USCRN <- function(start_date, end_date, dest_dir, DoNotOverwrite = TRUE) {
  # This function downloads precipitation data from the Boulder USCRN weather 
  # station at C1. It returns a list of the files it tried to download. By
  # default it will not download files that are already in the destination directory.
  # Arguments:
  #   start_date = the start date of tvan data in character form (or other form 
  #                that lubridate can coerce with its `year()` function)
  #   end_date = the end date of tvan data in character form (or other form 
  #              that lubridate can coerce with its `year()` function)
  #   dest_dir = the destination directory where the files will be downloaded
  #   DoNotOverwrite = should existing files with the same name be overwritten? If 
  #                 TRUE, files will not be overwritten, if FALSE, files will be 
                    #overwritten.
  
  
  require(lubridate)
  require(RCurl)
  
  # To do: replace this warning with a check for the tvan data
  message("Please note, end_date of USCRN data must not be less than the end_date of the tvan data.")
  
  
  # make dest_dir if it doesn't exist
  made_dir <- ifelse(!dir.exists(file.path(dest_dir)), 
                     dir.create(file.path(dest_dir), recursive = TRUE), FALSE)
  if (!made_dir) {
    writeLines("Data download directory not created, it already exists.")
  }
  
  # Create a list of urls - one for each year of data
  url_list <- vector(mode = "list", 
                     length = lubridate::year(end_date) - lubridate::year(start_date) + 1)
  
  file_list <- vector(mode = "list", 
                     length = lubridate::year(end_date) - lubridate::year(start_date) + 1)
  
  # get the names for each year (including unfinished partial years at the end)
  names(url_list) <- lubridate::year(seq(from = lubridate::ymd(as.Date(start_date)), 
                                         length.out = (lubridate::year(end_date) - 
                                                         lubridate::year(start_date) + 1),
                                         by = "years"))
  names(file_list) <- lubridate::year(seq(from = lubridate::ymd(as.Date(start_date)), 
                                         length.out = (lubridate::year(end_date) - 
                                                         lubridate::year(start_date) + 1),
                                         by = "years"))
  
  for (i in seq_along(url_list)) {
    url_list[[i]] <- paste0("https://www1.ncdc.noaa.gov/pub/data/uscrn/products/subhourly01/", names(url_list[i]), "/CRNS0101-05-", names(url_list[i]),"-CO_Boulder_14_W.txt")
  }
  
  
  # Check if url exists and if it does, download file
  for (i in seq_along(url_list)) {
  
    writeLines(paste0("Checking if ", url_list[[i]], " exists..."))
    if (!url.exists(url_list[[i]])) {
      stop(paste0("Url ", x, " is not accessible."))
    } else {
      writeLines("TRUE")
    }
    
    # Check if destination file already exists
      dest_fp <- paste0(dest_dir, "/CRNS0101-05-", 
                                 names(url_list[i]),"-CO_Boulder_14_W.txt")
      file_list[[i]] <- dest_fp
      
      if (file.exists(dest_fp) & DoNotOverwrite == TRUE) {
        writeLines(paste0(dest_fp, " already exits, skipping..."))
      } else { # if file doesn't exist or if overwrite is TRUE, download
        try(download.file(url = url_list[[i]], 
                      destfile = dest_fp))
      }
  }
  return(file_list)
}

# Function for reading in USCRN precip text files
read_USCRN_precip_data <- function(USCRN_precip_fp) {
  # This function reads in USCRN precipitation data files. It adds column
  # names and then it 1) collapses the time from 5-minute increments to half-
  # hourly by summing the precipitation over each 1/2-hour period; 2) Changes -9999
  # to NAs; and 3) selects only the local date, local time, and precpitation variables
  # for the final data frame. It returns the resulting dataframe.
  # Arguments:
  #     USCRN_precip_fp = file path to the USCRN text file you want to load
  
  
  # USCRN Fields and information can be found here:
  # https://www1.ncdc.noaa.gov/pub/data/uscrn/products/subhourly01/README.txt
  # Field#  Name                           Units
  # ---------------------------------------------
  # 1    WBANNO                         XXXXX
  # 2    UTC_DATE                       YYYYMMDD
  # 3    UTC_TIME                       HHmm
  # 4    LST_DATE                       YYYYMMDD
  # 5    LST_TIME                       HHmm
  # 6    CRX_VN                         XXXXXX
  # 7    LONGITUDE                      Decimal_degrees
  # 8    LATITUDE                       Decimal_degrees
  # 9    AIR_TEMPERATURE                Celsius
  # 10   PRECIPITATION                  mm
  # 11   SOLAR_RADIATION                W/m^2
  # 12   SR_FLAG                        X
  # 13   SURFACE_TEMPERATURE            Celsius
  # 14   ST_TYPE                        X
  # 15   ST_FLAG                        X
  # 16   RELATIVE_HUMIDITY              %
  # 17   RH_FLAG                        X
  # 18   SOIL_MOISTURE_5                m^3/m^3
  # 19   SOIL_TEMPERATURE_5             Celsius
  # 20   WETNESS                        Ohms
  # 21   WET_FLAG                       X
  # 22   WIND_1_5                       m/s
  # 23   WIND_FLAG                      X
  #              
  # ----------------------- Begin Function -------------------- # 
  
  require(dplyr)
  
  # read in text file
  writeLines(paste0("Reading in ", USCRN_precip_fp))
  precip <- read.table(USCRN_precip_fp, sep = "", 
                       colClasses = c(rep("character", times = 6), 
                                      rep("numeric", times = 7),
                                      "character",
                                      rep("numeric", times = 9)))
  # Assign column names
  names(precip) <- c("WBANNO", "UTC_DATE", "UTC_TIME", "LST_DATE", "LST_TIME",
                     "CRX_VN", "LONGITUDE", "LATITUDE", "AIR_TEMPERATURE", 
                     "PRECIPITATION", "SOLAR_RADIATION", "SR_FLAG", 
                     "SURFACE_TEMPERATURE", "ST_TYPE", "ST_FLAG",
                     "RELATIVE_HUMIDITY", "RH_FLAG", "SOIL_MOISTURE_5", 
                     "SOIL_TEMPERATURE_5", "WETNESS", "WET_FLAG", "WIND_1_5", 
                     "WIND_FLAG")
  
  # Clean data frame
  precip <- precip %>% 
    # Split local time string and convert to decimal time
    dplyr::mutate(UTC_TIME = gsub("(..)(..)", "\\1:\\2:00", UTC_TIME),
           cleanTime_UTC =
             strsplit(UTC_TIME, ":") %>%
             sapply(function(x){
               x <- as.numeric(x)
               x[1] + x[2]/60 + x[3]/(60*60)
             }), 
           decimalTime_UTC = floor(cleanTime_UTC * 2)/2) %>% 
    dplyr::mutate(LST_TIME = gsub("(..)(..)", "\\1:\\2:00", LST_TIME),
                  cleanTime_LST =
                    strsplit(LST_TIME, ":") %>%
                    sapply(function(x){
                      x <- as.numeric(x)
                      x[1] + x[2]/60 + x[3]/(60*60)
                    }), 
                  decimalTime_LST = floor(cleanTime_LST * 2)/2) %>%
    # select only columns used for precipitation and time stamp
    dplyr::select(UTC_DATE, UTC_TIME, LST_DATE, LST_TIME, 
                  cleanTime_UTC, decimalTime_UTC, 
                  cleanTime_LST, decimalTime_LST, PRECIPITATION) %>% 
    # set NAs from -9999
    dplyr::mutate_all(list(~na_if(., -9999))) %>% 
    # sum all precip events in each 1/2 period
    dplyr::group_by(UTC_DATE, decimalTime_UTC) %>%
    dplyr::mutate(PRECIP_TOT = sum(PRECIPITATION)) %>%
    # remove extra time steps
    dplyr::select(-PRECIPITATION, -LST_TIME, -UTC_TIME, 
                  -cleanTime_UTC, -cleanTime_LST) %>%
    unique() %>%
    # create 1/2-hourly time stamps
    dplyr::mutate(UTC_DATE = as.Date(UTC_DATE, format = "%Y%m%d"),
                  timestamp_UTC = as.POSIXct(paste0(UTC_DATE," 00:00:00"), 
                                             tz = "UTC") + 3600*decimalTime_UTC) %>%
    dplyr::mutate(LST_DATE = as.Date(LST_DATE, format = "%Y%m%d"),
                  timestamp_LST = as.POSIXct(paste0(LST_DATE," 00:00:00"),
                                             tz = "MST") + 3600*decimalTime_LST)
  
  return(precip)
}

# Function for downloading radiation data from Ameriflux
download_amflx <- function(dest_dir, username, 
                         site = "US-NR1", DescriptionOfDataUse,
                         DoNotOverwrite = TRUE,
                         verbose = FALSE) {
  # This function downloads radiation data from the Ameriflux webiste 
  # It returns a list of the files it tried to download. By default it will 
  # not download files that are already in the destination directory.
  # Arguments:
  #   dest_dir  -------------- the destination directory where the files will be 
  #                            downloaded
  #   username  -------------- the Ameriflux username of the user - this function 
  #                            will fail without a valid username.
  #   site  ------------------ the Ameriflux site to get the data from; defaults to 
  #                            US-NR1
  #   DescriptionOfDataUse --- the description to provide to Ameriflux for the intended
  #                            use of the data. If not provided by the user, the 
  #                            description will read:
  #
  #                            These data will be used as atmospheric forcings 
  #                            to run a local point-simulation for the alpine 
  #                            tundra at the Niwot Ridge LTER site.
  #
  #   DoNotOverwrite --------- should existing files with the same name be overwritten? 
  #                            If TRUE, files will not be overwritten, if FALSE, files
  #                            will be overwritten.
  #   verbose ---------------- Should the communication with the website be verbose?
  #                            default is FALSE.
  
  require(httr)
  require(jsonlite)
  require(RCurl)
  
  # Testing
  # site <- "US-NR1"
  # username <- amf_usr
  # dest_dir <- "~/Downloads/lter_flux/rad2"
  
  writeLines("Connecting with Ameriflux endpoint...")
  
  # NOTE THIS ENDPOINT MAY CHANGE
  ameriflux_endpoint <- "https://ameriflux-data.lbl.gov/AmeriFlux/DataDownload.svc/datafileURLs"
  
  if (missing(DescriptionOfDataUse)) {
    DescriptionOfDataUse = "These data will be used as atmospheric forcings to run a local point-simulation for the alpine tundra at the Niwot Ridge LTER site."
  }
  
  # Construct Payload request for ameriflux endpoint
  Payload <- paste0('{',
              '"username":"', username, '",',
              '"siteList":["', site, '"],',
              '"intendedUse": "Research - Land model/Earth system model",',
              '"description": "', DescriptionOfDataUse, '"',
                  '}') 
  
  # Get download information from Ameriflux endpoint
  if (verbose) {
    tmp <- httr::POST(url = ameriflux_endpoint,
                      body = Payload, verbose(), content_type_json())
  } else {
    tmp <- httr::POST(url = ameriflux_endpoint,
                      body = Payload, content_type_json())
  }
  
  # Check that the connection was successful
  if (tmp$status_code < 200 | tmp$status_code > 299) {
    stop(paste0("Attempt to connect to the website was not successful.\n",
                   "This may be because Ameriflux has changed its endpoint url \n",
                   "and you may need to contact Ameriflux support for an updated \n",
                   "address, or it may be due to a mistake in the request payload \n",
                   "syntax. Please check that the Ameriflux endpoing url and the \n", 
                   "payload syntax are valid. \n\n",
                   "Current endpoint: ", ameriflux_endpoint, "\n", 
                   "Current payload: ", Payload))
  } else {
    writeLines("Connection to Ameriflux successful.")
  }
  
  # extract content from the response
  r <- content(tmp) 
  
  # Check if the content is successfully received
  if (class(r) == "raw" | length(r$dataURLsList) == 0) {
    stop(paste0("No data was received from Ameriflux. Please check that your ",
                "username is valid and that both it and the site name are ",
                "spelled correctly."))
  }
  
  # Extract list of ftp urls
  url_list <- unlist(lapply(1:length(r$dataURLsList), 
                            function(x){r$dataURLsList[[x]]$URL}))
  
  file_list <- vector(mode = "list", 
                      length = length(url_list))

  # Notify user of the data policy prior to download
  message(paste0("Thank you for using Ameriflux data. Please be aware of the data \n",
                 "policy. By downloading this data you are acknowledging that you \n",
                 "have read and agree to that policy. \n\n",
                 "The following is how you described how you intend to use the data.\n\n",
                 "\tIntended Use: Research - Land model/Earth system model \n",
                 "\tDescription: These data will be used as atmospheric forcings \n",
                 "\tto run a local point-simulation for the alpine tundra at the \n",
                 "\tNiwot Ridge LTER site)\n\n",
                 "By downloading the data, the data contributors have been informed \n",
                 "of your use. If you are planning an in-depth analysis that may \n",
                 "result in a publication, please contact the data contributors \n",
                 "directly so they have the opportunity to contribute substantially \n",
                 "and become a co-author. \n\n", 
                 "The contact email for this site is: ", 
                 unlist(r$manifest$emailForSitePIs), "\n\n",
                 "You should also acknowledge Ameriflux in your presentations and \n",
                 "publications. Details about how this should be done can be found \n", 
                 "on the Ameriflux website. \n\n", 
                 "The full policy along with details about how to properly cite the \n",
                 "data can found here: \n",
                 "https://ameriflux.lbl.gov/data/data-policy/"))
  
  
  # make dest_dir if it doesn't exist
  made_dir <- ifelse(!dir.exists(file.path(dest_dir)), 
                     dir.create(file.path(dest_dir), recursive = TRUE), FALSE)
  
  writeLines("Downloading data...")
  
  if(!made_dir) {
    writeLines("Data download directory not created, it already exists.")
  }
  
  # Check if downloaded files already exist and if not, download file
  for (i in seq_along(url_list)) {
    
    # Check if destination file already exists
    dest_fp <- paste0(dest_dir, "/", basename(url_list[[i]]))
    file_list[[i]] <- dest_fp
    
    if (file.exists(dest_fp) & DoNotOverwrite == TRUE) {
      writeLines(paste0(dest_fp, " already exits, skipping..."))
    } else { # if file doesn't exist or if overwrite is TRUE, download
      # try(download.file(url = url_list[[i]], 
      #                   destfile = dest_fp, 
      #                   method = "curl"))
      try(GET(url = url_list[[i]], 
              write_disk(dest_fp, overwrite=FALSE), progress(), verbose()))
    }
  }
  return(unlist(file_list))
}

##############################################################################
# Read in L1 flux tower data product
##############################################################################
# Read in East & West tower
if (tower == "East" | tower == "Both") {
  # East data
  tvan_east <- read.table(file = east_data_fp, sep = "\t",
                          skip =  2, header = FALSE)
  tvan_east_names <- read.table(file = east_data_fp, sep = "\t",
                                header = TRUE, nrows = 1)
  tvan_east_units <- as.character(unname(unlist(tvan_east_names[1,])))
  
  colnames(tvan_east) <- names(tvan_east_names)
}
if (tower == "West" | tower == "Both") {
  # West data
  tvan_west <- read.csv(file = west_data_fp, sep = "\t",
                        skip =  2, header = FALSE)
  tvan_west_names <- read.table(file = west_data_fp, sep = "\t",
                                header = TRUE, nrows = 1)
  tvan_west_units <- as.character(unname(unlist(tvan_west_names[1,])))
  
  colnames(tvan_west) <- names(tvan_west_names)
}


# Get the start and end dates of the tvan data. If tower = "Both", 
# combine East and West data into one dataframe for convenience
if (tower == "Both") {
  tvan_east$Tower <- "East"
  tvan_west$Tower <- "West"
  
  tvan_all <- bind_rows(tvan_east, tvan_west) %>%
    mutate_all(list(~na_if(., -9999))) %>%
    mutate(date = as.Date(DoY - 1, origin = paste0(Year, "-01-01")),
           timestamp = as.POSIXct(paste0(date," 00:00:00"),
                                  format = "%Y-%m-%d %H:%M:%OS", 
                                  tz = "MST") + 3600*Hour) %>%
    group_by(Tower, Year, DoY) %>%
    mutate_at(vars(NEE:Ustar), list(daily_mean = mean), na.rm = TRUE) %>%
    select(date, timestamp, Year, DoY, Hour, Tower, everything())
  
  # Set a start/end date for the precip and radiation data based on the tvan data
  # make sure it's a round number or rEddyProc will complain
  start_date <- ceiling_date(min(tvan_all$timestamp, na.rm = TRUE), unit = "day")
  end_date <- floor_date(max(tvan_all$timestamp, na.rm = TRUE), unit = "day")
} else if (tower == "East") {
  tvan_east$Tower <- "East"
  # Set a start/end date for the precip and radiation data based on the tvan data
  start_date <- min(tvan_east$timestamp, na.rm = TRUE)
  end_date <- max(tvan_east$timestamp, na.rm = TRUE)
} else if (tower == "West") {
  tvan_west$Tower <- "West"
  # Set a start/end date for the precip and radiation data based on the tvan data
  start_date <- min(tvan_west$timestamp, na.rm = TRUE)
  end_date <- max(tvan_west$timestamp, na.rm = TRUE)
}


# Create a timeseries dataframe with the timestamps (this is in MST since start_date
# and end_date are in MST):
posix_complete <- as.data.frame(seq.POSIXt(start_date, end_date, by = "30 mins"))
colnames(posix_complete) <- "timestamp"
# get rid of first timestep, which is at midnight and not 00:30:00; it makes rEddyProc complain
posix_complete <- data.frame(timestamp = posix_complete[-1,])

##############################################################################
# Download Precipitation
##############################################################################
# Download precip data
# From here: https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nwt.416.10
writeLines("Downloading Saddle Precip data from EDI...")
saddle_precip_data_fp <- download_EDI(edi_id = saddle_precip_data, 
                                   dest_dir = paste0(DirDnld, "/precip_data"),
                                   getNewData = getNewData)


writeLines("Downloading C1 precipitation data from USCRN...")
USCRN_precip_data_fp <- download_USCRN(start_date = start_date, 
                                       end_date = end_date, 
                                       dest_dir = paste0(DirDnld, "/precip_data"),
                                       DoNotOverwrite = TRUE)


##############################################################################
# Handling Precip data
##############################################################################
# Saddle precip data must be corrected for blowing snow events, and extended to 
# half-hourly precip using Will's formula (see below for details).

writeLines("Reading in Saddle data...")
# Read in Saddle and USCRN Precip data; also collapse USCRN data into one dataframe
saddle_precip <- read.csv(saddle_precip_data_fp,
                          sep = ",", quot = '"', check.names = TRUE)

writeLines("Reading in C1 precipitation data from USCRN. This may take a while.")
USCRN_precip_list <- lapply(USCRN_precip_data_fp, read_USCRN_precip_data)
USCRN_precip <- plyr::rbind.fill(USCRN_precip_list) %>% 
  unique() # make sure to remove duplicates caused by aggregating to 30-minute time steps

# Check for duplicated time stamps - should be 0 (aka no TRUEs)
if (sum(duplicated(USCRN_precip$timestamp_UTC)) > 0) {
  warning("USCRN precipitation data still contains ",
          sum(duplicated(USCRN_precip$timestamp_UTC)),
          " duplicates!")
} else {
  writeLines(paste0("USCRN precipitation data has been loaded. ",
                    sum(duplicated(USCRN_precip$timestamp_UTC)), 
                    " duplicated timestamps have been detected."))
}


# Filter the precip data by exact start and end dates
saddle_precip <- saddle_precip %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= floor_date(start_date, unit = "day") & 
           date <= ceiling_date(end_date, unit = "day"))

USCRN_precip <- USCRN_precip %>%
  rename(date = LST_DATE) %>%
  mutate(timestamp_LST = as.POSIXct(timestamp_LST, tz = "MST")) %>%
  filter(timestamp_LST >= floor_date(start_date, unit = "day") &
           timestamp_LST <= ceiling_date(end_date, unit = "day"))

# Apply blowing snow correction to months of Oct-May Saddle data
# Due to blowing snow events where the belfort gauge has an oversampling of precipitation,
# it is recommended to add a correction for the precipitation total in the months Oct-May.
# The recommended correction for these events should be (0.39 * the recorded total). More
# information on this can be found in: 
# Williams, M.W., Bardsley, T., Rikkers, M., (1998) Overestimation of snow depth and inorganic nitrogen wetfall using NADP data, Niwot Ridge, Colorado. Atmospheric Environment 32 (22) :3827-3833 
writeLines("Applying blowing snow correction to Saddle precip data.")
saddle_precip <- saddle_precip %>%
  mutate(month = month(date),
         ppt_tot_corr = ifelse(month %in% c(10, 11, 12, 1, 2, 3, 4, 5), 
                               ppt_tot * 0.39, ppt_tot))

# Change any Nas or NaNs to zero
saddle_precip <- saddle_precip %>% 
  mutate(ppt_tot_corr = ifelse(is.na(ppt_tot_corr), 0, ppt_tot_corr))

USCRN_precip <- USCRN_precip %>% 
  mutate(PRECIP_TOT = ifelse(is.na(PRECIP_TOT), 0, PRECIP_TOT))


# Apply Will's algorithm for Precip data from paper:
# Use half-hourly precipitation recordfrom the U.S. Climate Reference Network (USCRN; data from https://www1.ncdc.noaa.gov/pub/data/uscrn/products/subhourly01/;), measured nearby (4 km) at the lower elevation(3050 m asl) C-1 site. Proportioanlly allocate the daily saddle precip measurements to the half-hourly precip record from USCRN. On days when Saddle record reports measurable precip, but the USCRN does not, distribute the daily saddle precip evenly across the day for model simulations. 
# Code modified from his TVAN_daily_ppt.R script
writeLines(paste0("Applying Will Wieder's algorithm for allocating daily Saddle ", 
                  "precipitation totals into 30-minute increments."))
Tvan_ppt <- saddle_precip$ppt_tot_corr
CRNS_ppt <- USCRN_precip$PRECIP_TOT
CRNS_date <- USCRN_precip$date
CRNS_mo <- month(USCRN_precip$date)
CRNS_hour <- USCRN_precip$decimalTime
CRNS_d    <- tapply(CRNS_ppt,  CRNS_date, sum)	# daily precip totals
CRNS_day  <- tapply(CRNS_date, CRNS_date, mean)	# num of days since 1970-01-01 - see date.mean()
CRNS_month <- tapply(CRNS_mo, CRNS_date, mean) # months 	
#------------------------------------------------------
# distribute Tvan ppt when observed in half-hourly CRNS
#------------------------------------------------------
ndays        <- length(Tvan_ppt)
nsteps       <- length(CRNS_ppt)
Tvan_fine    <- rep(NA, nsteps)
Tvan_note    <- rep(NA, nsteps)
Tvan_flag    <- rep(NA, ndays)
Tvan_flag_mo <- rep(NA, ndays)
Tvan_date    <- USCRN_precip$date # MST date
Tvan_hour    <- USCRN_precip$decimalTime_LST # MST hour
start        <- 1	

# code below does the following:
  # (0) if no daily precip at Tvan, add zeros to half hourly results 
  # (1) if precip at Tvan, but not recorded @ CRNS, distribute evenly in day and add 1 the flag
  # (2) if both precip at Tvan and CRNS, distribute Tvan in same proportion as CRNS 

for (d in 1:ndays) {
  end <- start + 47
  
  if (Tvan_ppt[d] == 0) {
    Tvan_fine[start:end] <- 0
    Tvan_note[start:end] <- 0
  } else if (CRNS_d[d] == 0){
    Tvan_fine[start:end] <- Tvan_ppt[d] / 48		
    Tvan_note[start:end] <- 1
    Tvan_flag[d]         <- 1
    Tvan_flag_mo[d]      <- CRNS_month[d]
  } else {
    temp_frac <- CRNS_ppt[start:end] / CRNS_d[d]
    Tvan_fine[start:end] <- Tvan_ppt[d] * temp_frac
    Tvan_note[start:end] <- 2
  }

  if (round(sum(Tvan_fine[start:end], na.rm = TRUE), digits = 7) != 
      round(sum(Tvan_ppt[d], na.rm = TRUE), digits = 7)) {
     warning(paste0("Running precip totals don't match at day ", d))
   }
  
  
  start <- end + 1
}

# Check that the total precip that fell at the saddle is the same as the total precip
# when allocated over 30-minute time steps
if (sum(Tvan_fine, na.rm=T) == sum(Tvan_ppt)) {
  writeLines(paste0("Total precip that fell at the Saddle (", sum(Tvan_ppt), 
                    ") matches the amount of total precip that has been ",
                    "allocated to the for the tvan data (", sum(Tvan_fine, na.rm=T), ")."))
} else {
  warning(paste0("Total precip that fell at the Saddle (", sum(Tvan_ppt), 
                 ") does NOT match the amount of total precip that has been ",
                 "allocated to the for the tvan data (", sum(Tvan_fine, na.rm=T), ")!"))
}



writeLines(paste0("Number of total days = ",ndays, " [", ddays(ndays), "]"))
writeLines(paste0("Number of days w/ precip at Tvan = ", 
                 length(Tvan_ppt[Tvan_ppt > 0]))) 
writeLines(paste0("Number of days with Tvan precip but w/o recorded CRNS precip = ", 
                  sum(Tvan_flag, na.rm = T)))
hist(Tvan_flag_mo, xlim = c(1,12),
     main = paste0("Montly frequency of days with Tvan precip but ",
                   "w/o recorded CRNS precip"),
     xlab = "Months"
     )

# Convert precip from mm/30 minutes into mm/s
Precip = Tvan_fine[1:nsteps] # mm every 30 minutes
PRECTmms <- Precip / (30*60) # mm/s

# Combine date and 1/2-hourly precip into one dataframe and add a timestamp
hlf_hr_precip <- data.frame(PRECTmms = PRECTmms, # mm/s
                            MST_HOUR = Tvan_hour[1:nsteps], # decimal hours
                            MST_DATE = Tvan_date[1:nsteps]) %>% # date
  mutate(timestamp = as.POSIXct(paste0(MST_DATE," 00:00:00"), tz = "MST") + 
           3600*MST_HOUR) %>%
  # fix date so that "0" hour readings are converted into 24
  mutate(MST_DATE = if_else(MST_HOUR == 0, MST_DATE - 1, MST_DATE), 
         MST_HOUR = if_else(MST_HOUR == 0.0, 24, MST_HOUR))

##############################################################################
# Download Radiation data
##############################################################################
writeLines("Downloading Ameriflux radiation data...")
rad_data_fp <- download_amflx(dest_dir = paste0(DirDnld, "/rad_data"),
                            username = amf_usr, verbose = TRUE)

# Check if the files have already been unzipped, if not, unzip the zip file
for (i in seq_along(rad_data_fp)) {
  if (grepl(".zip", basename(rad_data_fp[i]))) {
    writeLines(paste0("Unzipping ", rad_data_fp[i]))
    # check if the unzipped files exist
    unzip_list <- unzip(zipfile = rad_data_fp[i], 
                        exdir = dirname(rad_data_fp[i]),
                        overwrite = FALSE)
  }
}

amf_data_fp <- list.files(dirname(rad_data_fp[i]), 
                          full.names = TRUE,
                          pattern = "*.csv")

##############################################################################
# Handle Radiation data
##############################################################################
# Note: Radiation data comes from the Ameriflux NR-1 site. Currently this 
# data cannot be downloaded automatically and has to be downloaded by hand from
# the Ameriflux site after getting a user account: https://ameriflux.lbl.gov/data/download-data/

# For CLM we will pull out incoming shortwave (necessary) and incoming longwave (optional).
# The net radation is provided by the Tvan tower datasets. 
# The possible Ameriflux variables are:
# NETRAD_1_1_2       (W m-2): Net radiation (no QA/QC or gapfilling)
# NETRAD_PI_F_1_1_2   (W m-2): Net radiation (gapfilled by tower team)
# SW_IN_1_1_1       (W m-2): Shortwave radiation, incoming (no QA/QC or gapfilling)
# LW_IN_1_1_1        (W m-2): Longwave radiation, incoming (no QA/QC or gapfilling)
# SW_IN_PI_F_1_1_1       (W m-2): Shortwave radiation, incoming (gapfilled by tower team)
# LW_IN_PI_F_1_1_1        (W m-2): Longwave radiation, incoming (gapfilled by tower team)
# SW_OUT_1_1_1       (W m-2): Shortwave radiation, outgoing (no QA/QC or gapfilling)
# LW_OUT_1_1_1        (W m-2): Longwave radiation, outgoing (no QA/QC or gapfilling)
# SW_OUT_PI_F_1_1_1       (W m-2): Shortwave radiation, outgoing (gapfilled by tower team)
# LW_OUT_PI_F_1_1_1        (W m-2): Longwave radiation, outgoing (gapfilled by tower team)


writeLines("Reading in Ameriflux radiation data...")
# Load in Radiation data:
amf_data <- read.csv(file = amf_data_fp[2], 
                     skip = 2, 
                     header = TRUE, 
                     na.strings = "-9999",
                     as.is = TRUE)

# Select timestamps, and radiation variables
rad_data <- amf_data[,c("TIMESTAMP_START", "TIMESTAMP_END", 
                        "SW_IN_1_1_1", # also sometimes called Rg
                        "LW_IN_1_1_1", # also sometimes called FLDS
                        "SW_IN_PI_F_1_1_1", # also sometimes called Rg
                        "LW_IN_PI_F_1_1_1", # also sometimes called FLDS
                        "SW_OUT_1_1_1",
                        "LW_OUT_1_1_1",
                        "SW_OUT_PI_F_1_1_1",
                        "LW_OUT_PI_F_1_1_1",
                        "NETRAD_1_1_2",
                        "NETRAD_PI_F_1_1_2")]

rad_data$TIMESTAMP_START <- as.POSIXct(as.character(rad_data$TIMESTAMP_START), format = "%Y%m%d%H%M%OS", tz = "MST")
rad_data$TIMESTAMP_END <- as.POSIXct(as.character(rad_data$TIMESTAMP_END), format = "%Y%m%d%H%M%OS", tz = "MST")

# Subset the radiation data to the Tvan time period, reformat the times to get hours
# and dates, finally, select only the radiation, hour, and date variables. 
hlf_hr_rad <- rad_data %>%
  mutate(date = lubridate::date(TIMESTAMP_END)) %>% 
  filter(date >= floor_date(start_date, unit = "day") & 
           date <= floor_date(end_date, unit = "day")) %>%
  # Take reading from end of period, keep the date at midnight as the day before
  # to be consistent with other variables
  mutate(MST_HOUR = lubridate::hour(TIMESTAMP_END) + 
           lubridate::minute(TIMESTAMP_END)/60,
         MST_DATE = lubridate::date(TIMESTAMP_START)) %>% 
  # fix date so that "0" hour readings are converted into 24
  mutate(MST_HOUR = if_else(MST_HOUR == 0.0, 24, MST_HOUR)) %>%
  # Calculate net radiation from in/out radiation
  mutate(radNet = (SW_IN_PI_F_1_1_1 - SW_OUT_PI_F_1_1_1) + 
           (LW_IN_PI_F_1_1_1 - LW_OUT_PI_F_1_1_1)) %>%
  rename(Rg_usnr1 = SW_IN_PI_F_1_1_1, FLDS = LW_IN_PI_F_1_1_1,
         SW_OUT = SW_OUT_PI_F_1_1_1, LW_OUT = LW_OUT_PI_F_1_1_1,
         timestamp = TIMESTAMP_END) %>%
  select(timestamp, MST_DATE, MST_HOUR, Rg_usnr1, FLDS, radNet)

##############################################################################
# Combine flux and met data
##############################################################################
if (tower == "East" | tower == "Both") {
  # East tower
  tvan_east_tms <- tvan_east %>% 
    mutate_all(list(~na_if(., -9999))) %>%
    mutate(date = as.Date(DoY - 1, origin = paste0(Year, "-01-01")),
           timestamp = as.POSIXct(paste0(date," 00:00:00"), 
                                  format = "%Y-%m-%d %H:%M:%OS",
                                  tz = "MST") + 3600*Hour) 
}

if (tower == "West" | tower == "Both") {
  # West tower
  tvan_west_tms <- tvan_west %>% 
    mutate_all(list(~na_if(., -9999))) %>%
    mutate(date = as.Date(DoY - 1, origin = paste0(Year, "-01-01")),
           timestamp = as.POSIXct(paste0(date," 00:00:00"),
                                  format = "%Y-%m-%d %H:%M:%OS",
                                  tz = "MST") + 3600*Hour) 
}

# Join the flux data to the posix_complete date sequence
if (tower == "Both") {
  tmp_east <- left_join(posix_complete, tvan_east_tms, by = "timestamp") %>%
    mutate(Tower = "East")
  tmp_west <- left_join(posix_complete, tvan_west_tms, by = "timestamp") %>%
    mutate(Tower = "West")
  
  tvan_comb_tms <- bind_rows(tmp_east, tmp_west)
  tvan_tms <- tvan_comb_tms %>%
    # Fill in the DoY, Hour, Date, and Year that are NAs
    mutate(date = lubridate::date(timestamp)) %>% 
    # Take reading from end of period, keep the date at midnight as the day before
    # to be consistent with other variables
    mutate(Hour = lubridate::hour(timestamp) + 
             lubridate::minute(timestamp)/60,
           date = lubridate::date(timestamp)) %>% 
    # fix date so that "0" hour readings are converted into 24
    mutate(Hour = if_else(Hour == 0.0, 24, Hour),
           date = if_else(Hour == 24, date-1, date),
           DoY = yday(date), 
           Year = year(date))
  
} else if (tower == "West") {
  tmp_west <- left_join(posix_complete, tvan_west_tms, by = "timestamp") %>%
    mutate(Tower = "West")
  tvan_tms <- tmp_west %>%
    # Fill in the DoY, Hour, Date, and Year that are NAs
    mutate(date = lubridate::date(timestamp)) %>% 
    # Take reading from end of period, keep the date at midnight as the day before
    # to be consistent with other variables
    mutate(Hour = lubridate::hour(timestamp) + 
             lubridate::minute(timestamp)/60,
           date = lubridate::date(timestamp)) %>% 
    # fix date so that "0" hour readings are converted into 24
    mutate(Hour = if_else(Hour == 0.0, 24, Hour),
           date = if_else(Hour == 24, date-1, date),
           DoY = yday(date), 
           Year = year(date))
} else {
  tmp_east <- left_join(posix_complete, tvan_east_tms, by = "timestamp") %>%
    mutate(Tower = "East")
  tvan_tms <- tmp_east %>%
    # Fill in the DoY, Hour, Date, and Year that are NAs
    mutate(date = lubridate::date(timestamp)) %>% 
    # Take reading from end of period, keep the date at midnight as the day before
    # to be consistent with other variables
    mutate(Hour = lubridate::hour(timestamp) + 
             lubridate::minute(timestamp)/60,
           date = lubridate::date(timestamp)) %>% 
    # fix date so that "0" hour readings are converted into 24
    mutate(Hour = if_else(Hour == 0.0, 24, Hour),
           date = if_else(Hour == 24, date-1, date),
           DoY = yday(date), 
           Year = year(date))
}

writeLines("Combining precipitation, radiation, and Tvan data.")
# Combine dataframes by date and time
dataDf <- tvan_tms %>%
  left_join(hlf_hr_precip, by = c("Hour" = "MST_HOUR", "date" = "MST_DATE", 
                                  "timestamp" = "timestamp")) %>%
  left_join(hlf_hr_rad, by = c("Hour" = "MST_HOUR", "date" = "MST_DATE", 
                               "timestamp" = "timestamp")) %>%
  select(timestamp, date, Year, DoY, Hour, Tower, everything())

# Renaming of variables:
# FLDS - incident longwave (FLDS) (W/m^2)
# FSDS - incident shortwave (FSDS, or Rg) (W/m^2) # Check that these are the same as SW_IN/LW_IN
# PRECTmms - precipitation (PRECTmms = PRECTmms) (mm/s)
# PSRF - pressure at the lowest atmospheric level (PSRF = P) (kPa)
# RH - relative humidity at lowest atm level (RH = rH) (%)
# TBOT - temperature at lowest atm level (TBOT = Tair) (K)
# WIND - wind at lowest atm level (WIND = U) (m/s)
# NEE - net ecosystem exchange (NEE = NEE) (umolm-2s-1) 
# FSH - sensible heat flux (FSH = H) (Wm-2)
# EFLX_LH_TOT - latent heat flux (EFLX_LH_TOT = LE) (Wm-2)
# GPP - gross primary productivity (GPP) (umolm-2s-1)
# Rnet - net radiation (Rnet = Rn) (W/m^2)

##############################################################################
# Plot the un-gapfilled data
##############################################################################
if (makeplots == TRUE) {
  # needs ggplot and dplyr/tidyr
  # change data to longform
  # Necessary for model:
  # tbot, wind, rh, PSRF, FLDS, FSDS, PRECTmms
  
  getgaplength <- function(gap, y = "notgap") {
    res <- rle(gap == y)
    res_vec <- rep(res$values*res$lengths,res$lengths)
    return(res_vec)
  }
  
  # Find the minimum and maximum time stamps at which all required forcing variables have values
  min_gap_days <- 1 # how many days does a gap have to be at minimum to be plotted
  dataClm.forc.gaps <- dataDf %>%
    rename(TIMESTAMP = timestamp, EFLX_LH_TOT = LE, FSH = H, TBOT = Tair, RH = rH, 
                  WIND = U, PSRF = P, FSDS = Rg_usnr1) %>% 
    mutate_at(vars(TBOT, WIND, RH, PSRF, FLDS, FSDS, PRECTmms), list(gap = is.na)) %>%
    mutate(gap = TBOT_gap | WIND_gap | RH_gap | PSRF_gap | FLDS_gap | FSDS_gap | 
             PRECTmms_gap) %>%
    group_by(Tower) %>%
    mutate(gap = ifelse(gap == FALSE, "notgap", "gap"),
           ncontiguousgaps = getgaplength(gap, "gap")) %>%
    filter(gap == "gap") %>%
    select(TIMESTAMP, gap, ncontiguousgaps, Tower) %>%
    mutate(ndays = ncontiguousgaps/48,
           ncontiguousgaps = as.factor(ncontiguousgaps)) %>%
    group_by(Tower, ndays) %>%
    summarize(min = min(TIMESTAMP, na.rm = TRUE), 
              max = max(TIMESTAMP, na.rm = TRUE)) %>%
    arrange(desc(ndays)) %>%
    mutate(ndays = as.factor(round(ndays, digits = 2))) %>%
    mutate(yr1 = year(min),
           yr2 = year(max)) %>%
    rowwise() %>%
    mutate(years = paste0(seq(yr1, yr2), collapse = " | ")) %>%
    select(-yr1, -yr2)
  
  # Plot the required forcing variables
  dataClm.forc.plot <- dataDf %>%
    rename(TIMESTAMP = timestamp, EFLX_LH_TOT = LE, FSH = H, TBOT = Tair, RH = rH, 
           WIND = U, PSRF = P, FSDS = Rg_usnr1) %>%
    tidyr::pivot_longer(cols = !matches(c("TIMESTAMP", "date", "Year", "DoY", "Hour",
                                          "Tower")), 
                        names_to = "variable", 
                        values_to = "value") %>%
    filter(variable %in% c("TBOT", "WIND", "RH", "PSRF", "FLDS", "FSDS", "PRECTmms"))
  
  
  
  plot_gaps <- function(forcings, gaps, 
                        filteryears = NA, 
                        tower = NA, 
                        min_gap_days = 1,
                        highlightgaps = FALSE,
                        verbose = FALSE) {
    # if filteryear and tower are NA all years and both towers are plotted.
    # filteryear takes values of  either NA or a vector of character strings 
    # of years to plot
    # if highlightgaps == TRUE, gaps will be highlighted on plot
    # min_gaps_days is the minimum length in days of gaps to highlight
    forcings.plot <- forcings
    gaps.plot <- gaps %>%
      filter(as.numeric(as.character(ndays)) >= min_gap_days)
    if (nrow(gaps.plot) < 1) {highlightgaps = FALSE}

    title <- paste0("Gap-plots\n",
                    "Both towers \n", 
                    "Years: all")
    
    #### Filter forcing and gap datasets based on settings ####
    # create a custom title
    if (any(!is.na(filteryears)) & !is.na(tower)) { # filter towers and years
      forcings.plot <- forcings %>% 
        filter(Year %in% filteryears) %>%
        filter(Tower == tower) %>%
        # the following variables are the same in both towers
        filter(!(variable %in% c("FLDS", "FSDS", "PRECTmms")))
      gaps.plot <- gaps.plot %>%
        filter(grepl(paste0(filteryears, collapse = "|"), years)) %>%
        filter(Tower == tower)
      if (nrow(gaps.plot) < 1) {highlightgaps = FALSE}
      
      title <- paste0("Gap-plots\n",
                      "Tower: ", tower, "\n", 
                      "Years: ", paste0(filteryears, collapse = ", "))
      
    } else if (any(!is.na(filteryears))) { # filter only by years
      forcings.plot <- forcings %>% 
        filter(Year %in% filteryears)
      gaps.plot <- gaps.plot %>%
        filter(grepl(paste0(filteryears, collapse = "|"), years))
      if (nrow(gaps.plot) < 1) {highlightgaps = FALSE}
      
      title <- paste0("Gap-plots\n",
                      "Both towers \n",
                      "Years: ", paste0(filteryears, collapse = ", "))
      
    } else if (!is.na(tower)) { # filter only by tower
      forcings.plot <- forcings %>%
        filter(Tower == tower) %>%
        # the following variables are the same in both towers
        filter(!(variable %in% c("FLDS", "FSDS", "PRECTmms"))) 
      
      gaps.plot <- gaps.plot %>%
        filter(Tower == tower)
      if (nrow(gaps.plot) < 1) {highlightgaps = FALSE}
      
      title <- paste0("Gap-plots\n",
                      "Tower: ", tower, "\n", 
                      "Years: all")
    }
    
    # Tell the user what's happening
    writeLines(paste0("Plotting from ", min(forcings.plot$Year, na.rm = TRUE), " to ",
                      max(forcings.plot$Year, na.rm = TRUE)))
    if (verbose) {
      if (!is.na(tower)) {
        writeLines(paste0("Tower is ", tower))
      } else {
        writeLines(paste0("Plotting both towers"))
      }
      
      if (highlightgaps) {
        writeLines("Gaps will be highlighted")
        writeLines("Note: if a gap exeeds the boundary year, the x-axis will be", 
                   "modified so the entire gap is shown but points for that period ",
                   "will not be plotted.")
      } else { 
        writeLines("Gaps will not be highlighted")
      }
    }
    
    #### Plot the data ####
    forcing_gaps.plot <- ggplot(forcings.plot) +
      geom_point(aes(x = TIMESTAMP, y = value, color = Tower), alpha = 0.05) +
      facet_wrap(~variable, scales = "free_y", ncol = 1) +
      scale_color_discrete(name = "Tower") +
      guides(color = guide_legend(override.aes = list(alpha = 1),
                                  title.position = "top")) +
      theme(legend.position = "bottom") +
      ggtitle(title)
    
    # Highlight gaps on graphs
    if (highlightgaps) {
      forcing_gaps.plot <- forcing_gaps.plot +
        geom_rect(data = gaps.plot,
                  aes(xmin = min,
                      xmax = max,
                      ymin = -Inf, 
                      ymax = Inf,
                      fill = Tower), alpha = 0.3) +
        geom_vline(aes(xintercept = min), data = gaps.plot) +
        geom_vline(aes(xintercept = max), data = gaps.plot) +
        theme(legend.position = "bottom") +
        scale_fill_discrete(name = paste0("Gaps >", min_gap_days, " days")) +
        guides(fill = guide_legend(title.position = "top"))
    }
    return(forcing_gaps.plot)
  }
  
  plot_years <- c(min(dataClm.forc.plot$Year, na.rm = TRUE):max(dataClm.forc.plot$Year, na.rm = TRUE))
  plot_years <- set_names(plot_years)
  plot_years <- map(plot_years,
                    ~plot_gaps(forcings = dataClm.forc.plot, 
                               gaps = dataClm.forc.gaps, 
                               highlightgaps = TRUE, 
                               filteryears = .x, 
                               tower = NA,
                               min_gap_days = 7))
  
  plot_all_years <- plot_gaps(forcings = dataClm.forc.plot, 
                              gaps = dataClm.forc.gaps, 
                              highlightgaps = FALSE, 
                              filteryears = NA, 
                              tower = NA,
                              min_gap_days = 7)
  
  writeLines("Saving plots - this may take a while...")
  iwalk(plot_years, ~{
    suppressWarnings(
      ggsave(plot = .x, 
             filename = paste0(plots_dir,"/","yearly_gap_plots_",
                             .y, '.png'), 
             width = 10, 
             height = 5*7, 
             dpi = 150)
    )
    })
  
  forc.plot.out.name <- paste0(plots_dir,"/","all_years_gap_plots.png")
  
  ggsave(plot = plot_all_years, 
         filename = forc.plot.out.name,
         width = 10, 
         height = 5*7, 
         dpi = 150)
  
}

plots_dir

##############################################################################
# Gap-fill West tower with East tower 
##############################################################################
if (tower == "Both") {
  writeLines(paste0("Gap-filling ", basetower," tower data with data from the",
                    " other tower"))
  dataDf.wide <- dataDf %>% 
    select(all_of(c("timestamp", "date", "Year", "DoY", "Hour", "Tower", 
                     "NEE", "LE", "H", "Ustar", "Tair", "VPD", "rH", "U",
                     "P", "Rg_usnr1", "PRECTmms", "FLDS", "radNet", "Tsoil"))) %>%
    rename(Rg = Rg_usnr1) %>%
    mutate(BaseTower = ifelse(Tower == basetower, "base", "fill")) %>%
    # select(TIMESTAMP, date, Year, DoY, Hour, Tower, EFLX_LH_TOT, FSH,
    #                 TBOT, RH, WIND, PSRF, FSDS, FLDS, PRECTmms) %>%
    # for choice
    select(-Tower) %>%
    pivot_wider(names_from = BaseTower,
                values_from = c("NEE", "LE", "H", "Ustar", "Tair", "VPD",
                                "rH", "U", "P", "Rg", "PRECTmms",
                                "FLDS", "radNet", "Tsoil")) %>%
    # pivot_wider(names_from = Tower, 
    #           values_from = c("NEE", "LE", "H", "Ustar", "Tair", "VPD", 
    #                           "rH", "U", "P", "Rg", "PRECTmms",
    #                           "FLDS", "radNet", "Tsoil")) %>%
    select(!ends_with("_NA"))
  
  writeLines("Checking to make sure that tower timesteps line up correctly.")
  # convert posix_complete to UTC; then remove leap days
  #posix_complete$timestamp <- with_tz(posix_complete$timestamp, "UTC")
  # posix_complete_noleap <- posix_complete$timestamp[!grepl(".{4}-02-29", posix_complete$timestamp)]
  if (any(!(posix_complete$timestamp == dataDf.wide$timestamp))) {
    warning(paste0("At least one timestamp value is missing or out of bounds."))
  } else {
    writeLines(paste0("Timestamps are all present and line up correctly ",
                      "between \ntowers.",
                      "\nThere are ", nrow(dataDf.wide), 
                      " timestamps in total which is \n",
                      ddays(nrow(dataDf.wide)/48)))
  }
  
  #### Gap-fill "base" tower with "fill" tower data ####
  # we will create a flag variable to show which values were substituted
  # s = base tower was gapfilled with fill tower data
  # m = missing in both tower datasets
  # n = not missing; original west tower value was used
  gap_filled_from_twr <- dataDf.wide %>%
    mutate(
      # LH (Latent heat flux)
      LE = ifelse(is.na(LE_base), 
                           LE_fill, 
                           LE_base),
      LE_flag = ifelse(is.na(LE_base) & is.na(LE_fill), 
                                "m", 
                         ifelse(is.na(LE_base) & !is.na(LE_fill),
                                "s", "n")),
     # H (sensible heat flux)
     H = ifelse(is.na(H_base), 
                          H_fill, 
                          H_base),
     H_flag = ifelse(is.na(H_base) & is.na(H_fill), 
                               "m", 
                               ifelse(is.na(H_base) & !is.na(H_fill),
                                      "s", "n")),
     # Air Temperature (TBOT)
     Tair = ifelse(is.na(Tair_base), 
                          Tair_fill, 
                          Tair_base),
     Tair_flag = ifelse(is.na(Tair_base) & is.na(Tair_fill), 
                               "m", 
                               ifelse(is.na(Tair_base) & !is.na(Tair_fill),
                                      "s", "n")),
     # Relative humidity (rH)
     rH = ifelse(is.na(rH_base), 
                   rH_fill, 
                   rH_base),
     rH_flag = ifelse(is.na(rH_base) & is.na(rH_fill), 
                        "m", 
                        ifelse(is.na(rH_base) & !is.na(rH_fill),
                               "s", "n")),
     # Wind speed (U)
     U = ifelse(is.na(U_base), 
                   U_fill, 
                   U_base),
     U_flag = ifelse(is.na(U_base) & is.na(U_fill), 
                        "m", 
                        ifelse(is.na(U_base) & !is.na(U_fill),
                               "s", "n")),
     # Atmospheric pressure (P)
     P = ifelse(is.na(P_base), 
                   P_fill, 
                   P_base),
     P_flag = ifelse(is.na(P_base) & is.na(P_fill), 
                        "m", 
                        ifelse(is.na(P_base) & !is.na(P_fill),
                               "s", "n")),
     # Incident shortwave radiation (Rg_usnr1)
     Rg = ifelse(is.na(Rg_base), 
                   Rg_fill, 
                   Rg_base),
     Rg_flag = ifelse(is.na(Rg_base) & is.na(Rg_fill), 
                        "m", 
                        ifelse(is.na(Rg_base) & !is.na(Rg_fill),
                               "s", "n")),
     # Incident longwave radiation (FLDS) <- CHECK WITH WILL ON THIS ONE
     FLDS = ifelse(is.na(FLDS_base), 
                   FLDS_fill, 
                   FLDS_base),
     FLDS_flag = ifelse(is.na(FLDS_base) & is.na(FLDS_fill), 
                        "m", 
                        ifelse(is.na(FLDS_base) & !is.na(FLDS_fill),
                               "s", "n")),
     # Precipitation (PRECTmms)
     PRECTmms = ifelse(is.na(PRECTmms_base), 
                   PRECTmms_fill, 
                   PRECTmms_base),
     PRECTmms_flag = ifelse(is.na(PRECTmms_base) & is.na(PRECTmms_fill), 
                        "m", 
                        ifelse(is.na(PRECTmms_base) & !is.na(PRECTmms_fill),
                               "s", "n")),
     # Net Ecosystem Excahange (NEE)
     NEE = ifelse(is.na(NEE_base), 
                       NEE_fill, 
                       NEE_base),
     NEE_flag = ifelse(is.na(NEE_base) & is.na(NEE_fill), 
                            "m", 
                            ifelse(is.na(NEE_base) & !is.na(NEE_fill),
                                   "s", "n")), 
     # Ustar friction velocity (Ustar)
     Ustar = ifelse(is.na(Ustar_base), 
                  Ustar_fill, 
                  Ustar_base),
     Ustar_flag = ifelse(is.na(Ustar_base) & is.na(Ustar_fill), 
                       "m", 
                       ifelse(is.na(Ustar_base) & !is.na(Ustar_fill),
                              "s", "n")),
     # Net radiation (radNet)
     radNet = ifelse(is.na(radNet_base), 
                  radNet_fill, 
                  radNet_base),
     radNet_flag = ifelse(is.na(radNet_base) & is.na(radNet_fill), 
                       "m", 
                       ifelse(is.na(radNet_base) & !is.na(radNet_fill),
                              "s", "n")),
     # Soil Temperature (Tsoil)
     Tsoil = ifelse(is.na(Tsoil_base), 
                     Tsoil_fill, 
                     Tsoil_base),
     Tsoil_flag = ifelse(is.na(Tsoil_base) & is.na(Tsoil_fill), 
                          "m", 
                          ifelse(is.na(Tsoil_base) & !is.na(Tsoil_fill),
                                 "s", "n"))
     
     )
  #### Save Gap-filled outputs ####
  writeLines("Tower gap-filling complete. Saving data with flags...")
  dataDf <- gap_filled_from_twr %>% 
    select(!ends_with(c("base", "fill", "flag")))
  dataDf_flag <- gap_filled_from_twr %>% 
    select(!ends_with(c("base", "fill")))

  twr <- ifelse(tower == "Both", "both_towers", paste0(tower, "_tower"))
  flagged_fp <- paste0(DirOut, "/", "tvan_forcing_data_flagged_", 
                       twr, '_',lubridate::date(start_date),
                       '_',lubridate::date(end_date),".txt")
  write(paste0("# Flags: \n", 
               "# Base tower is: ", basetower, "\n",
               "# s = base tower was gapfilled with fill tower data \n",
               "# m = missing in both tower datasets \n",
               "# n = not missing; original west tower value was used"), 
        flagged_fp)
  suppressWarnings(
    write.table(dataDf_flag, flagged_fp, 
                sep = "\t", row.names = FALSE, 
                append = TRUE)
    )
  writeLines(paste0("Flagged data can be found here: ", flagged_fp))
}

##############################################################################
# Prepare file for ReddyProc
##############################################################################
# Change NA to -9999
dataDf[is.na(dataDf)] <- -9999

# #Convert time to ReddyProc format
# dataDf$Year <- lubridate::year(dataDf$TIMESTAMP) 
# dataDf$DoY <- lubridate::yday(dataDf$TIMESTAMP) 
# dataDf$Hour <- lubridate::hour(dataDf$TIMESTAMP) + lubridate::minute(dataDf$TIMESTAMP)/60
# 
# Remove timestamp and date
dataDf$timestamp <- NULL
dataDf$date <- NULL

# FLDS - incident longwave (FLDS) (W/m^2)
# FSDS - incident shortwave (FSDS) (W/m^2) # Check that these are the same as SW_IN/LW_IN
# PRECTmms - precipitation (PRECTmms = PRECTmms) (mm/s)
# PSRF - pressure at the lowest atmospheric level (PSRF = P) (kPa) - CONVERT TO kPa
# RH - relative humidity at lowest atm level (RH = rH) (%)
# TBOT - temperature at lowest atm level (TBOT = Tair) (K)
# WIND - wind at lowest atm level (WIND = U) (m/s)
# NEE - net ecosystem exchange (NEE = NEE) (umolm-2s-1) 
# FSH - sensible heat flux (FSH = H) (Wm-2)
# EFLX_LH_TOT - latent heat flux (EFLX_LH_TOT = LE) (Wm-2)
# GPP - gross primary productivity (GPP) (umolm-2s-1)
# Rnet - net radiation (Rnet = Rn/Rg) (W/m^2)
# Ustar - friction velocity
# Tsoil


#Vector of units for each variable
unitDf <- c("Year" = "--", "DoY" = "--", "Hour" = "--", "LE" = "Wm-2", "H" = "Wm-2",
            "Tair" = "degC", 
            "rH" = "%", "U" = "ms-1", "P" = "kPa",  "Rg" = "Wm-2", "FLDS" = "Wm-2", 
            "PRECTmms" = "mms-1", 
            "NEE" = "umolm-2s-1", "Ustar" = "ms-1", "radNet" = "Wm-2", 
            "Tsoil" = "degC")

#Set the output data column order based off of the units vector
dataDf <- data.table::setcolorder(dataDf, names(unitDf))

#Create filename
twr <- ifelse(tower == "Both", "both_towers", paste0(tower, "_tower"))

fileOut <- paste0(DirOut,"/","tvan_forcing_data_", 
                  twr, '_',lubridate::date(start_date),
                  '_',lubridate::date(end_date),'.txt')


h1 <- paste(names(unitDf), collapse = "\t")
h2 <- paste(unitDf, collapse = "\t")

#Output data in ReddyProc format
conFile <- file(fileOut, "w")
#write the variable names header
writeLines(text = c(h1,h2), sep = "\n", con = conFile)
#write the variable units header
#writeLines(text = unitDf, sep = "\t", con = conFile)
#Write output in tab delimited format
write.table(x = dataDf, file = conFile, sep = "\t", row.names = FALSE, col.names = FALSE)

#Close file connection
close(conFile)

##############################################################################
# ReddyProc Gap-filling workflow
##############################################################################
EddyData.F <- fLoadTXTIntoDataframe(fileOut)

#Threshold bounds to prevent rH > 100%
EddyData.F$rH[EddyData.F$rH > 100] <- 100
#Threshold bounds to prevent Rg (FSDS) < 0
EddyData.F$Rg[EddyData.F$Rg < 0] <- 0
EddyData.F$Rg[EddyData.F$Rg   > 1200 ] <- 1200
#Threshold bounds to prevent NEE > 100
EddyData.F$NEE[EddyData.F$NEE > 100] <- NA
#Threshold bounds to prevent NEE < -100
EddyData.F$NEE[EddyData.F$NEE < -100] <- NA

#+++ If not provided, calculate VPD from TBOT and RH
EddyData.F <- cbind(EddyData.F,VPD = fCalcVPDfromRHandTair(EddyData.F$rH, 
                                                           EddyData.F$Tair))

#+++ Add time stamp in POSIX time format
EddyDataWithPosix.F <- fConvertTimeToPosix(EddyData.F, 'YDH', Year = 'Year', 
                                           Day = 'DoY', Hour = 'Hour', tz = "MST")


#+++ Initalize R5 reference class sEddyProc for processing of eddy data
#+++ with all variables needed for processing later
EddyProc.C <- sEddyProc$new(twr, EddyDataWithPosix.F, 
                            c('NEE','Rg','Tair','VPD','rH','LE','H','Ustar','P', 
                              'FLDS','U', 'PRECTmms', 'radNet', 'Tsoil'))
#Set location information
EddyProc.C$sSetLocationInfo(LatDeg = latSite, LongDeg = lonSite, TimeZoneHour = -6)

#+++ Fill gaps in variables with MDS gap filling algorithm (without prior ustar filtering)
# Note, this also takes a long time to complete!
EddyProc.C$sMDSGapFill('NEE', FillAll = TRUE) #Fill all values to estimate flux uncertainties
EddyProc.C$sMDSGapFill('LE', FillAll = TRUE)
EddyProc.C$sMDSGapFill('H', FillAll = TRUE)
EddyProc.C$sMDSGapFill('Ustar', FillAll = TRUE)
EddyProc.C$sMDSGapFill('Tair', FillAll = FALSE)  
EddyProc.C$sMDSGapFill('VPD',    FillAll = FALSE) 
EddyProc.C$sMDSGapFill('rH',     FillAll = FALSE) 
EddyProc.C$sMDSGapFill('U', FillAll = FALSE) # wind
EddyProc.C$sMDSGapFill('PRECTmms', FillAll = FALSE) 
EddyProc.C$sMDSGapFill('P', FillAll = FALSE) 
EddyProc.C$sMDSGapFill('FLDS', FillAll = FALSE) 
EddyProc.C$sMDSGapFill('Rg', FillAll = FALSE) 
EddyProc.C$sMDSGapFill('radNet', FillAll = FALSE) 
EddyProc.C$sMDSGapFill('Tsoil', FillAll = FALSE) 
EddyProc.C$sMRFluxPartition()

#+++ Export gap filled and partitioned data to standard data frame
FilledEddyData.F <- EddyProc.C$sExportResults()

#Grab just the filled data products
dataClm <- FilledEddyData.F[,grep(pattern = "_f$", x = names(FilledEddyData.F))]

#Grab the POSIX timestamp
dataClm$DateTime <- EddyDataWithPosix.F$DateTime - lubridate::minutes(30) # putting back to original position

names(dataClm) <- gsub("_f", "", names(dataClm))

#Convert degC to K for temperature
dataClm$Tair <- dataClm$Tair + 273.15
attributes(obj = dataClm$Tair)$units <- "K"
#Convert kPa to Pa for pressure
dataClm$P <- dataClm$P * 1000.0
attributes(obj = dataClm$P)$units <- "Pa"

#Create tower height measurement field
dataClm$ZBOT <- rep(2,nrow(dataClm))

#Year month combination for data filtering
dataClm$yearMon <- paste0(year(dataClm$DateTime), "-", 
                          sprintf("%02d", month(dataClm$DateTime)))

##############################################################################
# Plotting and identifying gaps left in data after gapfilling
##############################################################################

if (makeplots == TRUE) {
  # needs ggplot and dplyr/tidyr
  # change data to longform
  # Necessary for model:
  # tbot, wind, rh, PSRF, FLDS, FSDS, PRECTmms
  
  getgaplength <- function(gap, y = "notgap") {
    res <- rle(gap == y)
    res_vec <- rep(res$values*res$lengths,res$lengths)
    return(res_vec)
  }
  
  # Find the minimum and maximum time stamps at which all required forcing variables have values
  dataClm.forc.gaps <- dataClm %>%
    rename(EFLX_LH_TOT = LE, FSH = H, TBOT = Tair, RH = rH, 
           WIND = U, PSRF = P, FSDS = Rg) %>%
    mutate_at(vars(TBOT, WIND, RH, PSRF, FLDS, FSDS, PRECTmms), list(gap = is.na)) %>%
    mutate(gap = TBOT_gap | WIND_gap | RH_gap | PSRF_gap | FLDS_gap | FSDS_gap | 
             PRECTmms_gap) %>%
    mutate(gap = ifelse(gap == FALSE, "notgap", "gap"),
           ncontiguousgaps = getgaplength(gap, "gap")) %>%
    filter(gap == "gap") %>%
    select(DateTime, gap, ncontiguousgaps) %>%
    mutate(ndays = ncontiguousgaps/48,
           #ndays = as.factor(ndays),
           ncontiguousgaps = as.factor(ncontiguousgaps)) %>%
    group_by(ndays) %>%
    summarize(min = min(DateTime, na.rm = TRUE), 
              max = max(DateTime, na.rm = TRUE)) %>%
    arrange(desc(ndays)) %>%
    mutate(ndays = as.factor(round(ndays, digits = 2))) %>%
    mutate(yr1 = year(min),
           yr2 = year(max)) %>%
    rowwise() %>%
    mutate(years = paste0(seq(yr1, yr2), collapse = " | ")) %>%
    select(-yr1, -yr2)
  
  
  # Plot the required forcing variables
  dataClm.forc.plot <- dataClm %>%
    rename(EFLX_LH_TOT = LE, FSH = H, TBOT = Tair, RH = rH, 
           WIND = U, PSRF = P, FSDS = Rg) %>%
    tidyr::pivot_longer(cols = !matches(c("DateTime", "yearMon")), 
                        names_to = "variable", 
                        values_to = "value") %>%
    filter(variable %in% c("TBOT", "WIND", "RH", "PSRF", "FLDS", "FSDS", "PRECTmms"))
  
  
  
  plot_gaps <- function(forcings, gaps, 
                        filteryears = NA,
                        min_gap_days = 1,
                        highlightgaps = FALSE,
                        verbose = FALSE) {
    # if filteryear and tower are NA all years and both towers are plotted
    # filteryear is either NA or a vector of character strings of years to plot
    # if highlightgaps == TRUE, gaps will be highlighted on plot
    # min_gaps_days is the minimum length in days of gaps to highlight
    forcings.plot <- forcings %>%
      mutate(Year = year(DateTime))
    gaps.plot <- gaps
    title <- paste0("Gap-plots for both towers and all years")
    if (any(!is.na(filteryears))) {
      forcings.plot <- forcings.plot %>% 
        filter(Year %in% filteryears)
      gaps.plot <- gaps.plot %>%
        filter(grepl(paste0(filteryears, collapse = "|"), years))
      
      title <- paste0("Gap-plots for gap-filled data: year(s) ",
                      paste0(filteryears, collapse = ", "))
      
    }
    
    writeLines(paste0("Plotting from ", min(forcings.plot$Year), " to ",
                      max(forcings.plot$Year)))
    
    if (verbose) {
      if (!is.na(tower)) {
        writeLines(paste0("Tower is ", tower))
      } else {
        writeLines(paste0("Plotting both towers"))
      }
      
      if (highlightgaps) {
        writeLines("Gaps will be highlighted")
        writeLines("Note: if a gap exeeds the boundary year, the x-axis will be", 
                   "modified so the entire gap is shown but points for that period ",
                   "will not be plotted.")
      } else { 
        writeLines("Gaps will not be highlighted")
      }
    }
    
    forcing_gaps.plot <- ggplot(forcings.plot) +
      geom_point(aes(x = DateTime, y = value), alpha = 0.05) +
      facet_wrap(~variable, scales = "free_y", ncol = 1) +
      ggtitle(title)
    
    if (nrow(gaps.plot) == 0) {
      highlightgaps <- FALSE
    }
    if (highlightgaps) {
      forcing_gaps.plot <- forcing_gaps.plot +
        geom_rect(data = gaps.plot,
                  aes(xmin = min,
                      xmax = max,
                      ymin = -Inf, 
                      ymax = Inf), alpha = 0.3) +
        geom_vline(aes(xintercept = min), data = gaps.plot) +
        geom_vline(aes(xintercept = max), data = gaps.plot) +
        theme(legend.position = "none") +
        scale_fill_discrete(name = paste0("Gaps >", min_gap_days, " days"))
    }
    return(forcing_gaps.plot)
  }
  
  plot_years <- c(min(year(dataClm.forc.plot$DateTime), 
                      na.rm = TRUE):max(year(dataClm.forc.plot$DateTime), 
                                        na.rm = TRUE))
  plot_years <- set_names(plot_years)
  plot_years <- map(plot_years,
                    ~plot_gaps(forcings = dataClm.forc.plot, 
                               gaps = dataClm.forc.gaps, 
                               highlightgaps = TRUE, 
                               filteryears = .x))
  
  iwalk(plot_years, ~{
    ggsave(plot = .x, 
           filename = paste0(plots_dir,"/",.y,
                             '_yearly_gap_plots_postgapfilling.png'), 
           width = 10, 
           height = 5*7, 
           dpi = 150)
  })
  
  
  plot_all_years <- plot_gaps(forcings = dataClm.forc.plot, 
                              gaps = dataClm.forc.gaps, 
                              highlightgaps = TRUE, 
                              filteryears = NA)
  
  forc.plot.out.name <- paste0(plots_dir,"/",
                               lubridate::date(dataClm$DateTime[1]),'_',
                               lubridate::date(tail(dataClm$DateTime, n = 1)),
                               '_required_forcing_postgapfilling.png')
  ggsave(plot = plot_all_years, 
         filename = forc.plot.out.name,
         width = 10, 
         height = 5*7, 
         dpi = 150)
  
}

##############################################################################
# Prepare 4 different precipitation regimes for the different vegetation communities
##############################################################################
# There are several vegetation communities at Niwot and they all see slightly 
# different precipitation regimes. (See Wieder et al. 2017). We will modify the 
# precipitation inputs based on Table 1 in Wieder et al. 2017

# | Community         | Snow (% relative to observations)       |
# | ----------------- | --------------------------------------  |
# | Fellfield (FF)    | 10, but 25 during March, April and May  |
# | Dry meadow (DM)   | 10, but 25 during March, April and May  |
# | Moist meadow (MM) | 100                                     |
# | Wet meadow (WM)   | 75 + runoff simulated from moist meadow |
# | Snowbed (SB)      | 200                                     |

dataClm_veg_communities <- dataClm %>% 
  mutate(month = month(DateTime),
         PRECTmms_FF = ifelse(Tair >= 273.15, PRECTmms,
                              ifelse(month %in% c(3,4,5), PRECTmms * 0.25,
                                     PRECTmms*0.1)),
         PRECTmms_DM = ifelse(Tair >= 273.15, PRECTmms,
                              ifelse(month %in% c(3,4,5), PRECTmms * 0.25,
                                     PRECTmms*0.1)),
         PRECTmms_MM = PRECTmms,
         PRECTmms_WM = ifelse(Tair >= 273.15, PRECTmms, PRECTmms*0.75),
         PRECTmms_SB = ifelse(Tair >= 273.15, PRECTmms, PRECTmms*2)) %>%
  select(-month)

# Add in simulated runoff from mm to wm:
if (simulated_runoff_present) {
  simulated_runoff <- read.csv(file = simulated_runoff_fp)
  colnames(simulated_runoff) <- names(simulated_runoff)
}

names(dataClm_veg_communities)
# convert runoff time to DateTime
simulated_runoff$time = as.POSIXct(simulated_runoff$time,tz='UTC')
simulated_runoff$time = round(simulated_runoff$time, 'min')

# add runoff to precipitation for wetmeadow
if(simulated_runoff_present){
  dataClm_veg_communities = dataClm_veg_communities %>%
    left_join(simulated_runoff, by = c("DateTime" = "time")) %>%
    mutate(PRECTmms_WM = PRECTmms_WM + QRUNOFF ) %>%
    select(-QRUNOFF)
}


# Write out modified precipitation data
twr <- ifelse(tower == "Both", "both_towers", paste0(tower, "_tower"))
precip_mods_fp <- paste0(DirOut, "/", "tvan_forcing_data_precip_mods_", 
                     twr, '_',lubridate::date(start_date),
                     '_',lubridate::date(end_date),".txt")

# ADD UNITS
dataClm_veg_communities_units <- c("NEE" = "umolm-2s-1", 
                                   "LE" = "Wm-2", 
                                   "H" =  "Wm-2", 
                                   "Ustar" = "ms-1",
                                   "Tair" = "K",
                                   "VPD" = "kPa",
                                   "rH" = "%",
                                   "U" = "ms-1",
                                   "PRECTmms" = "mms-1",
                                   "PRECTmms_FF" = "mms-1",
                                   "PRECTmms_DM" = "mms-1",
                                   "PRECTmms_MM" = "mms-1",
                                   "PRECTmms_WM" = "mms-1",
                                   "PRECTmms_SB" = "mms-1",
                                   "P" = "Pa",
                                   "FLDS" = "Wm-2",
                                   "Rg" = "Wm-2", 
                                   "radNet" = "Wm-2", 
                                   "Tsoil" = "degC",
                                   "GPP" = "umolm-2s-1",
                                   "DateTime" = "-",
                                   "yearMon" = "-", 
                                   "ZBOT" = "-")

# Reorder the units to match the order of dataClm_veg_communities
dataClm_veg_communities_units <- dataClm_veg_communities_units[names(dataClm_veg_communities)]
dataClm_veg_communities_units.df <- rbind(dataClm_veg_communities_units)
rownames(dataClm_veg_communities_units.df) <- NULL
  
write.table(dataClm_veg_communities_units.df, precip_mods_fp, 
            sep = "\t", row.names = FALSE)  
write.table(dataClm_veg_communities, precip_mods_fp, 
              sep = "\t", row.names = FALSE, append = TRUE, col.names = FALSE)  
##############################################################################
# Write output to CLM
##############################################################################
write_to_clm <- function(dataClm, veg_community = NA, verbose = FALSE) {
  # dataClm = the gap-filled data subsetted according to the precipitation
  #           regime you want
  # veg_community = one of "FF", "DM", "MM", "WM", or "SB" specifying the 
  #               vegetation community you want to simulate, if NA, original
  #               precip values are used
  #
  
  # Set up for vegetation choice
  veg_community_list <- c("fell_field", "dry_meadow", "moist_meadow", "wet_meadow",
                          "snow_bed")
  names(veg_community_list) <- c("FF", "DM", "MM", "WM","SB")
  if (is.na(veg_community)) { # original precip
    dataClm <- dataClm %>% 
      select(!ends_with(c("_FF", "_DM", "_MM", "_WM", "_SB")))
    vegcom <- "original"
  } else { # specific vegetation community
    precip_col_name <- paste0("PRECTmms_", veg_community)
    dataClm$PRECTmms <- dataClm[,precip_col_name]
    dataClm <- dataClm %>%
      select(!ends_with(c("_FF", "_DM", "_MM", "_WM", "_SB")))
    vegcom <- veg_community_list[veg_community]
  }
  
  
  #Define missing value fill
  mv <- -9999.  
  
  #Set of year/month combinations for netCDF output
  setYearMon <- unique(dataClm$yearMon)
  
  for (m in setYearMon) {
    #m <- setYearMon[10] #for testing
    Data.mon <- dataClm[dataClm$yearMon == m,]
    timeStep <- seq(0,nrow(Data.mon)-1,1)
    time     <- timeStep/48
    #endStep  <- startStep + nsteps[m]-1
    
    if (verbose) {
      print(paste(m,"Data date =",Data.mon$DateTime[1], "00:00:00"))
      names(Data.mon)
    }
    
    #NetCDF output filename
    fileOutNcdf <- paste(DirOut,"/",vegcom, "/",m,".nc", sep = "")
    if (verbose) {
      print(fileOutNcdf)
    }
    
    veg_com_dir <- paste0(DirOut,"/",vegcom)
    if(!dir.exists(veg_com_dir)) dir.create(veg_com_dir, recursive = TRUE)
    #sub(pattern = ".txt", replacement = ".nc", fileOut)
    
    
    
    # define the netcdf coordinate variables (name, units, type)
    lat  <- ncdf4::ncdim_def("lat","degrees_north", as.double(latSite), create_dimvar=TRUE)
    lon <- ncdf4::ncdim_def("lon","degrees_east", as.double(lonSite), create_dimvar=TRUE)
    
    #Variables to output to netCDF
    time <- ncdf4::ncdim_def("time", paste("days since",Data.mon$DateTime[1], "00:00:00"),
                             vals=as.double(time),unlim=FALSE, create_dimvar=TRUE,
                             calendar = "gregorian")
    LATIXY  <- ncdf4::ncvar_def("LATIXY", "degrees N", list(lat), mv,
                                longname="latitude", prec="double")
    LONGXY  <- ncdf4::ncvar_def("LONGXY", "degrees E", list(lon), mv,
                                longname="longitude", prec="double")
    FLDS  <- ncdf4::ncvar_def("FLDS", "W/m^2", list(lon,lat,time), mv,
                              longname="incident longwave (FLDS)", prec="double")
    FSDS  <- ncdf4::ncvar_def("FSDS", "W/m^2", list(lon,lat,time), mv,
                              longname="incident shortwave (FSDS)", prec="double")
    PRECTmms <- ncdf4::ncvar_def("PRECTmms", "mm/s", list(lon,lat,time), mv,
                                 longname="precipitation (PRECTmms)", prec="double")
    PSRF  <- ncdf4::ncvar_def("PSRF", "Pa", list(lon,lat,time), mv,
                              longname="pressure at the lowest atmospheric level (PSRF)", prec="double")
    RH    <- ncdf4::ncvar_def("RH", "%", list(lon,lat,time), mv,
                              longname="relative humidity at lowest atm level (RH)", prec="double")
    TBOT  <- ncdf4::ncvar_def("TBOT", "K", list(lon,lat,time), mv,
                              longname="temperature at lowest atm level (TBOT)", prec="double")
    WIND  <- ncdf4::ncvar_def("WIND", "m/s", list(lon,lat,time), mv,
                              longname="wind at lowest atm level (WIND)", prec="double")
    ZBOT  <- ncdf4::ncvar_def("ZBOT", "m", list(lon,lat,time), mv,
                              longname="observational height", prec="double")
    NEE <- ncdf4::ncvar_def("NEE", "umolm-2s-1", list(lon,lat,time), mv,
                              longname="net ecosystem exchange", prec="double")
    FSH  <- ncdf4::ncvar_def("FSH", "Wm-2", list(lon,lat,time), mv,
                             longname="sensible heat flux", prec="double")
    EFLX_LH_TOT  <- ncdf4::ncvar_def("EFLX_LH_TOT", "Wm-2", list(lon,lat,time), mv,
                                     longname="latent heat flux", prec="double")
    GPP <- ncdf4::ncvar_def("GPP", "umolm-2s-1", list(lon,lat,time), mv,
                            longname="gross primary productivity", prec="double")
    Rnet  <- ncdf4::ncvar_def("Rnet", "W/m^2", list(lon,lat,time), mv,
                              longname="net radiation", prec="double")
    
    #Create the output file
    ncnew <- ncdf4::nc_create(fileOutNcdf, list(LATIXY,LONGXY,FLDS,FSDS,PRECTmms,RH,PSRF,TBOT,WIND,ZBOT,FSH,EFLX_LH_TOT,NEE,GPP,Rnet))
    
    
    # Write some values to this variable on disk.
    ncdf4::ncvar_put(ncnew, LATIXY, latSite)
    ncdf4::ncvar_put(ncnew, LONGXY, lonSite)
    ncdf4::ncvar_put(ncnew, FLDS, Data.mon$FLDS)
    ncdf4::ncvar_put(ncnew, FSDS, Data.mon$Rg)
    ncdf4::ncvar_put(ncnew, RH,   Data.mon$rH)
    ncdf4::ncvar_put(ncnew, PRECTmms, Data.mon$PRECTmms)
    ncdf4::ncvar_put(ncnew, PSRF, Data.mon$P)
    ncdf4::ncvar_put(ncnew, TBOT, Data.mon$Tair)
    ncdf4::ncvar_put(ncnew, WIND, Data.mon$U)
    ncdf4::ncvar_put(ncnew, ZBOT, Data.mon$ZBOT)
    ncdf4::ncvar_put(ncnew, NEE, Data.mon$NEE)
    ncdf4::ncvar_put(ncnew, FSH, Data.mon$H)
    ncdf4::ncvar_put(ncnew, EFLX_LH_TOT, Data.mon$LE)
    ncdf4::ncvar_put(ncnew, GPP, Data.mon$GPP)
    ncdf4::ncvar_put(ncnew, Rnet, Data.mon$radNet)
    #add attributes
    # ncdf4::ncatt_put(ncnew, time,"calendar", "gregorian" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, FLDS,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, FSDS,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, RH  ,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, PRECTmms,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, PSRF,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, TBOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, WIND,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, ZBOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, NEE,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, FSH,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, EFLX_LH_TOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, GPP,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, Rnet,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, 0, "veg_community_type", veg_community_list[veg_community],prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, 0, "created_on",date()       ,prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, 0, "created_by","Hannah Holland-Moritz",prec=NA,verbose=FALSE,definemode=FALSE )
    ncdf4::ncatt_put(ncnew, 0, "created_from",fileOut        ,prec=NA,verbose=FALSE,definemode=FALSE )

    ncdf4::ncatt_put(ncnew, 0, "created_with", "flow.lter.clm.R",prec=NA,verbose=FALSE,definemode=FALSE )
    
    
    #Close Netcdf file connection
    ncdf4::nc_close(ncnew)
    #Add step
    #startStep <- endStep + 1
    #Remove not needed variables
    remove(time, timeStep, fileOutNcdf, ncnew, Data.mon,
           FLDS,FSDS,RH,PRECTmms,PSRF,TBOT,WIND,ZBOT)
  } #End of monthloop
}


# Prepare file for CLM simulations - convert to UTC and filter out leapdays
dataClm_veg_communities_modelready <- dataClm_veg_communities %>%
  # Convert time into UTC
  mutate(timestamp_UTC = with_tz(DateTime, tzone = "UTC"),
         date = as.Date(timestamp_UTC),
         Hour = lubridate::hour(timestamp_UTC) + 
           lubridate::minute(timestamp_UTC)/60) %>%
  # Remove leap years
  filter(!grepl(".{4}-02-29", date)) %>% 
  # Fix Hours, date, DoY, and Year; Hour is 0.5-24.0; Adjust date accordingly
  # get new doy now that leap years are filtered out
  mutate(Hour = if_else(Hour == 0.0, 24, Hour),
         date = if_else(Hour == 24, date - 1, date),
         Year = year(date),
         DoY = yday(date),
         DoY = ifelse(leap_year(Year) & (yday(date) > 59),
                      (yday(date) - 1), yday(date))) %>%
  # Remove MST timestamp and replace it with UTC timestamp; also remove other
  # extraneous time indicators
  select(-DateTime, -date, -Hour, -Year, -DoY) %>%
  rename(DateTime = timestamp_UTC) %>%
  # overwrite yearMon with updated timezone yearMon
  mutate(yearMon = paste0(year(DateTime), "-", 
                          sprintf("%02d", month(DateTime)))) 


# Create NC files
community_list <- c("Fell Field", "Dry Meadow", "Moist Meadow", "Wet Meadow",
                      "Snow Bed", "Original Precipitation")
names(community_list) <- c("FF", "DM", "MM", "WM","SB", NA)
for (i in seq_along(community_list)) {
  writeLines(paste0("Writing .nc files for ", 
                    community_list[i], "..."))
  write_to_clm(dataClm = dataClm_veg_communities_modelready, 
               veg_community = names(community_list[i]))
}

print(DirOut)
print('The met (.nc) forcings for Tvan are ready to be used! Time to run CLM')

