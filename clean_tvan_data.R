################################################################################
# @title T-Van supplemental cleaning

# @author
# Hannah Holland-Moritz \email{Hannah.HollandMoritz@colorado.edu}

# @description 
# Supplemental cleaning of tvan data in preparation for the tvan workflow that goes beyond
# the cleaning done in tvan_L1_preprocess.R script
# Notably, this script adds Saddle air temperature to Tvan Tair results when the sensors 
# seemed to go bad (2016-2018)

# changelog and author contributions / copyrights
# Hannah Holland-Moritz (2020-06-09)
#   Writing the initial script
# Will Wieder (2020-08-24)
#   Minor changes and suggestions 
##############################################################################

##############################################################################
# Dependencies
##############################################################################
rm(list = ls())

packReq <- c("dplyr", "ggplot2", "tidyr", "lubridate", "EML", "xts")

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
# Should plots be made?
makeplots <- TRUE # TRUE = default

#### Output Options ####
# The output directory for the script. It is recommended but not required that this
# be the same directory that holds the Reddyproc-ready files produced by 
# tvan_L1_preprocess.R
user = 'wwieder'
if (user ==  'wwieder') {
  DirOutBase <- paste0("~/Desktop/Working_files/Niwot/Tvan_out_new")
  east_data_fp <- paste0(DirOutBase,"/Reddy_proc_readyData/tvan_East_2007-08-29_09-00-00_to_2021-01-05_12-30-00_flux_P_reddyproc.txt")
  # The location of the west tvan data filepath, use "", if tower = "East"
  west_data_fp <- paste0(DirOutBase,"/Reddy_proc_readyData/tvan_West_2007-05-09_19-00-00_to_2021-01-04_00-30-00_flux_P_reddyproc.txt")
}  else {
  DirOutBase <- "~/Downloads/Tvan_out/Reddy_proc_readyData"
  east_data_fp <- paste0(DirOutBase,"/tvan_East_2007-08-29_09-00-00_to_2020-04-09_18-00-00_flux_P_reddyproc.txt")
  west_data_fp <-  paste0(DirOutBase,"/tvan_West_2007-05-09_19-00-00_to_2020-08-11_07-30-00_flux_P_reddyproc.txt")
}

#### Tower Use Options ####
# What tvan tower should be used?
tower <- "Both" # Options are "East", "West", or "Both"
# if "Both" the both towers will be processed at once


##############################################################################
# Static workflow parameters - these are unlikely to change
##############################################################################

#Append the site to the base output directory
DirOut <- paste0(DirOutBase, "/supp_filtering")
plots_dir <- paste0(DirOut, "/plots")

#Check if directory exists and create if not
if(!dir.exists(DirOut)) dir.create(DirOut, recursive = TRUE)
if(!dir.exists(plots_dir) & makeplots) dir.create(plots_dir, recursive = TRUE)

# the EDI id for hourly meterological data from the saddle weather stations
saddle_met_data <- "57" # NWT LTER EDI id


##############################################################################
# Helper functions - for downloading and loading data
##############################################################################
# Define a helpful filtering function 
filtering_function <- function(data, variable, start_date, end_date,
                               upper_cutoff, lower_cutoff) {
  # function for filtering points that can't be solved by one general cutoff
  # if a cutoff is NA, then no cutoff will be applied, otherise the cutoff will
  # be applied only within the starting and ending dates
  # 
  # Variable definitions:
  # data --------- the dataframe to be filtered
  # variable ----- the variable to be filtered
  # start_date --- the start date for the filtering period
  # end_date ----- the end date for the filtering period
  # upper_cutoff - the upper cutoff for filtering data; if NA, no cutoff will be applied
  # lower_cutoff - the lower cutoff for filtering data; if NA, no cutoff will be applied
  require(dplyr)

  # testing
  # print(names(data))
  # print(head(data[[variable]]))
  
  # if there is both an upper and lower cutoff
  if ( !is.na(upper_cutoff) & !is.na(lower_cutoff)) {
    
    # print("if 1")
    # filter the data
    data[[variable]] <- ifelse(as.Date(data[["timestamp"]]) > as.Date(start_date) & 
                                 as.Date(data[["timestamp"]]) < as.Date(end_date) & 
                                 (data[[variable]] > upper_cutoff |
                                    data[[variable]] < lower_cutoff), 
                               NaN, data[[variable]])
  } 
  # if there is no upper cutoff
  else if ( is.na(upper_cutoff)) {
    # print("if 2")
    # Filter the data
    data[[variable]] <- ifelse(as.Date(data[["timestamp"]]) > as.Date(start_date) & 
                                 as.Date(data[["timestamp"]]) < as.Date(end_date) & 
                                 data[[variable]] < lower_cutoff, 
                               NaN, data[[variable]])
  } 
  # if there is no lower_cutoff
  else if ( is.na(lower_cutoff)) {
    # print("if 3")
    # Filter the data
    data[[variable]] <- ifelse(as.Date(data[["timestamp"]]) > as.Date(start_date) & 
                                 as.Date(data[["timestamp"]]) < as.Date(end_date) & 
                                 data[[variable]] > upper_cutoff, 
                               NaN, data[[variable]])
  }
  
  
  return(list(variable = data[[variable]]))
  
}

plot_comparison_by_years <- function(year, data, outdir, var, plot_filetype='.pdf') {
  # Variable definitions
  # year ---------- The year of data to plot
  # data ---------- the dataframe containing the variable in both filtered and 
  #                 unfiltered states; Also needs a column "Filtered" to specify
  #                 which state the variable is in.
  # outdir -------- the location that the plots should be saved to
  # var ----------- the variable (character string format) being plotted
  # plot_filetype - the filtype (pdf or png) that should be saved, default to pdf
  
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  # Subset data by year and variable
  p.plot <- data %>% 
    dplyr::filter(!is.na(!!dplyr::sym(var))) %>% 
    dplyr::mutate(Year = lubridate::year(timestamp)) %>%
    filter(Year == year) %>%
    dplyr::select(timestamp, Year, all_of(var), Filtered, Tower)
  xmin <- as.POSIXct(paste0(year, "-01-01 00:00:00"), tz = "MST")
  xmax <- as.POSIXct(paste0(year + 1, "-01-01 00:00:00"), tz = "MST")
  # only plot years that have data
  if ( nrow(p.plot) > 2) {
    p <- ggplot(p.plot, aes_string(x = "timestamp", y = var)) +
      geom_point(aes_string(color = "Tower"), alpha = 0.3) + 
      facet_wrap(Year ~ Filtered, scales = "free") +
      ylab(paste0(var, " (", tvan_units[var], ")")) +
      xlim(xmin, xmax) +
      ggtitle(paste0(tower, " ", var, " (", tvan_units[var], ")", " - ", year))
    #p
    ggsave(paste0(outdir,"/", tower, "_", var, "_", year, plot_filetype),
           plot = p, width = 14, height = 7)
  }
}

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
  myeml<-xml2::read_xml(paste0("https://portal.edirepository.org/nis/metadataviewer?packageid=",

                                 packageid,
                                 "&contentType=application/xml")) %>% EML::read_eml()
}

# Function for downloading from EDI
download_EDI <- function(edi_id, dest_dir, getNewData = TRUE) {
  # This section heavily borrowed from Sarah Elmendorf's generic_timeseries_workflow.R script
  # https://github.com/NWTlter/long-term-trends/blob/master/plotting_scripts/generic_timeseries_workflow.R
  
  # Depends on getCurrentVersion() and getEML()
  library(dplyr)
  library(EML)
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
    
    # Download the files if they need to be downloaded; extract the column classes from the 
    # attributes and save the file paths and column classes of the data in the dest_dir 
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


# define a helpful plotting function
# # Helpful plotting
# data.frame(time = flux_P[["time"]], filtvariable = tmp$variable,
#            variable = flux_P[["tc200"]]) %>%
#   pivot_longer(contains("variable"), names_to = "OriginalFiltered",
#                values_to = "value") %>%
#   filter(as.Date(time) > as.Date("2017-10-01") &
#            as.Date(time) < as.Date("2017-12-31")) %>%
#   ggplot(aes(x = time, y = value, color = OriginalFiltered)) +
#   geom_point(alpha = 0.5)


##############################################################################
# Read in L1 flux tower data product
##############################################################################
# Read in East & West tower
if (tower == "East" | tower == "Both") {
  # East data
  tvan_east <- read.table(file = east_data_fp, sep = "\t",
                          skip =  2, header = FALSE)
  tvan_east_names <- read.table(file = east_data_fp, sep = "\t",
                                header = TRUE, nrows = 1, stringsAsFactors = FALSE)
  tvan_east_units <- as.character(unname(unlist(tvan_east_names[1,])))
  
  colnames(tvan_east) <- names(tvan_east_names)
  tvan_units <- tvan_east_names
}

if (tower == "West" | tower == "Both") {
  # West data
  tvan_west <- read.table(file = west_data_fp, sep = "\t",
                        skip =  2, header = FALSE)
  tvan_west_names <- read.table(file = west_data_fp, sep = "\t",
                                header = TRUE, nrows = 1, stringsAsFactors = FALSE)
  tvan_west_units <- as.character(unname(unlist(tvan_west_names[1,])))
  
  colnames(tvan_west) <- names(tvan_west_names)
  tvan_units <- tvan_west_names
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


# Create a timeseries dataframe with the timestamps:
posix_complete <- as.data.frame(seq.POSIXt(start_date, end_date, by = "30 mins"))
colnames(posix_complete) <- "timestamp"
# get rid of first timestep, which is at midnight and not 00:30:00; it makes rEddyProc complain
posix_complete <- data.frame(timestamp = posix_complete[-1,])

##############################################################################
# Organize data for plotting and filtering
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

# Convert rh to % 
if (tower == "East" | tower == "Both") {
  tvan_east_tms$rH <- tvan_east_tms$rH*100
}

if (tower == "West" | tower == "Both") {
  tvan_west_tms$rH <- tvan_west_tms$rH*100
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

# Save a version of the unfiltered data
tvan_tms.unfilt <- tvan_tms
tvan_tms.unfilt$Filtered <- "unfiltered"

##############################################################################
# Filter problem spots
##############################################################################
if (tower == "East" | tower == "Both") {
  
  #### Tair ####
  # Remove "flatlining" after september 2017-2019
  tmp <- filtering_function(data = tvan_east_tms, variable = "Tair",
                            start_date = "2017-07-01", end_date = "2019-03-30",
                            upper_cutoff = NA, lower_cutoff = 20)
  tvan_east_tms$Tair <- tmp$variable
  
  # data.frame(time = tvan_east_tms[["timestamp"]], filtvariable = tmp$variable,
  #            variable = tvan_east_tms[["Tair"]]) %>%
  #   pivot_longer(contains("variable"), names_to = "OriginalFiltered",
  #                values_to = "value") %>%
  #   filter(as.Date(time) > as.Date("2014-01-01") &
  #            as.Date(time) < as.Date("2020-12-31")) %>%
  #   ggplot(aes(x = time, y = value, color = OriginalFiltered)) +
  #   geom_point(alpha = 0.5)
  

  #### VPD ####
  # Remove "flatlining" after september 2017-2018
  tmp <- filtering_function(data = tvan_east_tms, variable = "VPD",
                            start_date = "2017-07-01", end_date = "2019-03-31",
                            upper_cutoff = NA, lower_cutoff = 20)
  
  tvan_east_tms$VPD <- tmp$variable
  
  # data.frame(time = tvan_east_tms[["timestamp"]], filtvariable = tmp$variable,
  #            variable = tvan_east_tms[["VPD"]]) %>%
  #   pivot_longer(contains("variable"), names_to = "OriginalFiltered",
  #                values_to = "value") %>%
  #   filter(as.Date(time) > as.Date("2017-01-01") &
  #            as.Date(time) < as.Date("2020-12-31")) %>%
  #   ggplot(aes(x = time, y = value, color = OriginalFiltered)) +
  #   geom_point(alpha = 0.5)
}

if (tower == "West" | tower == "Both") {
  
  #### Tair ####
  # Remove wildly low data at the end of 2011; only remove points that exceed the range
  # captured by the East tower, since the east and west tower have pretty good concordance.
  
  tmp <- filtering_function(data = tvan_west_tms, variable = "Tair",
                            start_date = "2011-11-28", end_date = "2012-01-01",
                            upper_cutoff = NA, lower_cutoff = -25)
  tvan_west_tms$Tair <- tmp$variable
  
  # data.frame(time = tvan_west_tms[["timestamp"]], filtvariable = tmp$variable,
  #            variable = tvan_west_tms[["Tair"]], 
  #            eastvariable = tvan_east_tms[match(tvan_west_tms$timestamp, 
  #                                               tvan_east_tms$timestamp), "Tair"]) %>%
  #   pivot_longer(contains("variable"), names_to = "OriginalFiltered",
  #                values_to = "value") %>%
  #   filter(as.Date(time) > as.Date("2011-01-01") &
  #            as.Date(time) < as.Date("2016-12-31")) %>%
  #   ggplot(aes(x = time, y = value, color = OriginalFiltered)) +
  #   geom_point(alpha = 0.5)
  
  # Remove "flatlining" after september 2016-2018
  tmp <- filtering_function(data = tvan_west_tms, variable = "Tair",
                            start_date = "2016-09-10", end_date = "2019-02-03",
                            upper_cutoff = NA, lower_cutoff = 20)
  tvan_west_tms$Tair <- tmp$variable
  
  # data.frame(time = tvan_west_tms[["timestamp"]], filtvariable = tmp$variable,
  #            variable = tvan_west_tms[["Tair"]]) %>%
  #   pivot_longer(contains("variable"), names_to = "OriginalFiltered",
  #                values_to = "value") %>%
  #   filter(as.Date(time) > as.Date("2014-01-01") &
  #            as.Date(time) < as.Date("2020-12-31")) %>%
  #   ggplot(aes(x = time, y = value, color = OriginalFiltered)) +
  #   geom_point(alpha = 0.5)
  
  
  #### VPD ####
  # Note: the data in 2011 don't look very different from the East tower, so I've
  # left them in even though Tair was used to calculate VPD here.
  # Remove "flatlining" after september 2016-2018
  tmp <- filtering_function(data = tvan_west_tms, variable = "VPD",
                            start_date = "2016-09-10", end_date = "2019-01-31",
                            upper_cutoff = NA, lower_cutoff = 20)
  tvan_west_tms$VPD <- tmp$variable
  
  # data.frame(time = tvan_west_tms[["timestamp"]], filtvariable = tmp$variable,
  #            variable = tvan_west_tms[["VPD"]]) %>%
  #   pivot_longer(contains("variable"), names_to = "OriginalFiltered",
  #                values_to = "value") %>%
  #   filter(as.Date(time) > as.Date("2018-01-01") &
  #            as.Date(time) < as.Date("2020-12-31")) %>%
  #   ggplot(aes(x = time, y = value, color = OriginalFiltered)) +
  #   geom_point(alpha = 0.5)
  
  #### rH ####
  # Remove wildly high data at the end of 2011 (probably caused by Tair)
  tmp <- filtering_function(data = tvan_west_tms, variable = "rH",
                            start_date = "2011-11-20", end_date = "2011-12-31",
                            upper_cutoff = 1.0, lower_cutoff = NA)
  tvan_west_tms$rH <- tmp$variable
  
  # data.frame(time = tvan_west_tms[["timestamp"]], filtvariable = tmp$variable,
  #            variable = tvan_west_tms[["Tair"]]) %>%
  #   pivot_longer(contains("variable"), names_to = "OriginalFiltered",
  #                values_to = "value") %>%
  #   filter(as.Date(time) > as.Date("2016-01-01") &
  #            as.Date(time) < as.Date("2019-12-31")) %>%
  #   ggplot(aes(x = time, y = value, color = OriginalFiltered)) +
  #   geom_point(alpha = 0.5)
}



##############################################################################
# Prepare data for plotting
##############################################################################
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

##############################################################################
# Plot data by years
##############################################################################
if (makeplots) {
  # Plot data
  tvan_tms.plot <- tvan_tms %>%
    mutate(Filtered = "filtered") %>%
    bind_rows(tvan_tms.unfilt) %>%
    mutate(Filtered = factor(Filtered, levels = c("unfiltered", "filtered"))) %>%
    select(timestamp, date, Year, DoY, Hour, Tower, Filtered, everything())
  
  # get sequence of years:
  year_list <- lubridate::year(seq.POSIXt(min(tvan_tms.plot$time), 
                                          max(tvan_tms.plot$time), by = "1 year"))
  
  # Plot by year for each variable
  for (i in 8:ncol(tvan_tms.plot)) {
    # Don't bother plotting any variable that is all NAs after filtering
    if (!all(is.na(tvan_tms.plot[[i]]))) {
      print(paste0("Plotting ", names(tvan_tms.plot)[i]))
      y_name = names(tvan_tms.plot)[i]
      yr_plots_dir <- paste0(plots_dir, "/",y_name, "_yearly_plots")
      ifelse(!dir.exists(file.path(yr_plots_dir)),
             dir.create(file.path(yr_plots_dir)), FALSE)
      writeLines(paste0("Saving yearly plots for ", y_name))
      lapply(year_list, plot_comparison_by_years,
             data = tvan_tms.plot,
             outdir = yr_plots_dir,
             var = y_name)
    }
    
  }
}

##############################################################################
# Read in Saddle Met data to replace Tair
##############################################################################
# Download Saddle Met data
message(paste0("Downloading Saddle Met data, please cite: \n",
               "Morse, J. and M. Losleben. 2019. Climate data for saddle data loggers (CR23X and CR1000), 2009 - ongoing, hourly. ver 3. Environmental Data Initiative. https://doi.org/10.6073/pasta/4f416341d978376c0205c86bc88d90ba (Accessed ",Sys.Date(), ")"))

saddle_met_data_fp <- download_EDI(edi_id = saddle_met_data, 
                                   dest_dir = paste0(DirOut,"/saddle_met_data"),
                                   getNewData = TRUE)

colclasses <- gsub("Date", "character", saddle_met_data_fp$colclasses)
sadd_met <- read.csv(file = saddle_met_data_fp$csv)#, 
#                     colClasses = colclasses, na.strings='NaN')


##############################################################################
# Filter saddle met data and compare with Tvan Tair
##############################################################################
sadd_met$date_time_start <- as.POSIXct(as.character(sadd_met$date_time_start),
                                       tz = "MST", format = "%Y-%m-%d %H:%M")
#sadd_met %>% 
#  filter(date_time_start > as.POSIXct("2019-01-01 00:00:00", tz = "MST")) %>%
#  select(date_time_start, airtemp_avg)

saddle_met_sub <- subset(sadd_met, date_time_start >= start_date)

# 2019 data should be mean of  airtemp_hmpX_avg, where  flags are OK

saddle_met_sub1 <- saddle_met_sub %>%
  filter(date_time_start > as.POSIXct("2019-01-01 00:00:00", tz = "MST")) %>%
  select(date_time_start, airtemp_avg, 
         airtemp_hmp1_avg, flag_airtemp_hmp1_avg,
         airtemp_hmp2_avg, flag_airtemp_hmp2_avg,
         airtemp_hmp3_avg, flag_airtemp_hmp3_avg) %>%
  mutate(airtemp_hmp1_avg = ifelse(flag_airtemp_hmp1_avg %in% c("m", "mq", "q"), 
                                   NA, airtemp_hmp1_avg)) %>%
  mutate(airtemp_hmp2_avg = ifelse(flag_airtemp_hmp2_avg %in% c("m", "mq", "q"), 
                                   NA, airtemp_hmp2_avg)) %>%
  mutate(airtemp_hmp3_avg = ifelse(flag_airtemp_hmp3_avg %in% c("m", "mq", "q"), 
                                   NA, airtemp_hmp3_avg)) %>%
  group_by(date_time_start) %>%
  summarise(airtemp_avg = mean(c(airtemp_hmp1_avg,airtemp_hmp2_avg,airtemp_hmp3_avg), na.rm=T))# %>%
  #select(date_time_start, airtemp_avg, airtemp_hmp1_avg)

plot(saddle_met_sub1)

# Remove missing or questionable observations and then remove the flag columns
saddle_met_sub2 <- saddle_met_sub %>%
  filter(date_time_start < as.POSIXct("2019-01-01 00:00:00", tz = "MST")) %>%

  select(date_time_start, airtemp_avg, flag_airtemp_avg) %>%
  mutate(airtemp_avg = ifelse(flag_airtemp_avg %in% c("m", "mq", "q"), 
                              NA, airtemp_avg)) %>%
  select(date_time_start, airtemp_avg)

#plot(saddle_met_sub2)
saddle_met_sub2 = bind_rows(saddle_met_sub2, saddle_met_sub1)
# Join airtemp_avg & airtemp_hmp_avg, usinig later for year > 2019
plot(saddle_met_sub2)

## Interploate met data from hourly to 1/2 hourly
sadd_met_alltime <- posix_complete %>%
  left_join(saddle_met_sub2, by = c("timestamp" = "date_time_start"))
met.xts <- xts(sadd_met_alltime, 
               order.by = sadd_met_alltime$timestamp, 
               tzone = "MST")
met.xts$airtemp_avg <- na.fill(met.xts$airtemp_avg, c(NA, "extend", NA))
met.df <- as.data.frame(met.xts) %>%
  mutate_all(as.character) %>%
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "MST"),
         airtemp_avg = as.numeric(airtemp_avg))

if (makeplots) {
  # Combine with Tvan for comparison plot
  # here the .pdf is GIANT, stick with .png
  sadd_met_compare.plot <- tvan_tms %>% 
    left_join(met.df, by = c("timestamp")) %>%
    select(timestamp, Year, DoY, Hour, Tair, airtemp_avg) %>%
    pivot_longer(matches(c("Tair", "airtemp_avg")),
                 names_to = "TemperatureSource", 
                 values_to = "AirTemperature") %>%
    mutate(TemperatureSource = ifelse(TemperatureSource == "Tair", "Tvan", "Saddle"))
  
  sadd_tvan_temp_comp <- ggplot(sadd_met_compare.plot,
                                aes(x = timestamp, y = AirTemperature)) +
    geom_point(aes(color = TemperatureSource), alpha = 0.03)
  
  ggsave(sadd_tvan_temp_comp, filename = paste0(DirOut, "/plots/sadd_tvan_temp_comp.png"),
         device = 'png', width = 10, height =  5, dpi = 150)
}

cortest.comp <- tvan_tms %>% 
  left_join(saddle_met_sub2, by = c("timestamp" = "date_time_start")) %>%
  select(timestamp, Tair, airtemp_avg)
cor.test(cortest.comp$Tair, cortest.comp$airtemp_avg)

# They seem pretty similar! We will fill gaps in the tvan tower data with the
# Saddle data
##############################################################################
# Use Saddle Air Temperature to fill gaps in Tvan data after 2016
##############################################################################
tvan_tms <- tvan_tms %>% 
  left_join(met.df, by = c("timestamp")) %>% 
  mutate(Tair = ifelse(is.na(Tair) & timestamp > "2016-01-01",
                       airtemp_avg, Tair)) %>%
  select(-airtemp_avg)

##############################################################################
# Output the data
##############################################################################
# set up dates for saving data
time1 <- gsub(":", "-", gsub(" ", "_", lubridate::with_tz(min(tvan_tms$time), tz = "MST")))
time2 <- gsub(":", "-", gsub(" ", "_", lubridate::with_tz(max(tvan_tms$time), tz = "MST")))
period <- paste0(time1, "_to_", time2)

if (tower == "East") {
  tvan_twr_east <- tvan_tms
}

if (tower == "West") {
  tvan_twr_west <- tvan_tms
}


if (tower == "Both") {
  # Filter by tower to get east and west
  tvan_twr_east <- tvan_tms %>% 
    filter(Tower == "East")
  
  tvan_twr_west <- tvan_tms %>%
    filter(Tower == "West") 
}

if (exists("tvan_twr_east")) {
  # Select columns and rename them according to ReddyProc rules
  # Note no radiation in this data!
  tvan_twr_east <- tvan_twr_east[,names(tvan_east_names)]
  
  colnames(tvan_twr_east) <- names(tvan_east_names)
  
  tvan_twr_east_units <- rbind(tvan_east_names, tvan_twr_east)
  
  # Write out tvan and met data
  # used `cleaned` to help clarify what was done to the reddyproc results
  writeLines(paste0("Saving ReddyProc-ready files to ",
                    paste0(DirOut, "/", 
                           "tvan_", "East","_", period, "_flux_P_reddyproc_cleaned.txt")))
  
  write.table(tvan_twr_east_units, 
              file = paste0(DirOut, "/", 
                            "tvan_", "East","_", period, "_flux_P_reddyproc_cleaned.txt"),
              row.names = FALSE, sep = "\t")
}

if (exists("tvan_twr_west")) {
  # Select columns and rename them according to ReddyProc rules
  # Note no radiation in this data!
  tvan_twr_west <- tvan_twr_west[,names(tvan_west_names)]
  
  colnames(tvan_twr_west) <- names(tvan_west_names)
  
  tvan_twr_west_units <- rbind(tvan_west_names, tvan_twr_west)
  
  # Write out tvan and met data
  writeLines(paste0("Saving ReddyProc-ready files to ",
                    paste0(DirOut, "/", 
                           "tvan_", "West","_", period, "_flux_P_reddyproc_cleaned.txt")))
  
  write.table(tvan_twr_west_units, 
              file = paste0(DirOut, "/", 
                            "tvan_", "West","_", period, "_flux_P_reddyproc_cleaned.txt"),
              row.names = FALSE, sep = "\t")
}

print('Finished cleaning Tvan Data, now create .nc files using prepare_forcings_for_clm.R')
