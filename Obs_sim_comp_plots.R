rm(list = ls())

##############################################################################
# Dependencies
##############################################################################

#Call the R HDF5 Library
packReq <- c("magrittr","EML", "dplyr", "ggplot2", 
             "purrr", "tidyr", "lubridate","RCurl", "cowplot")

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

# 1) Base directory for all files
# 2) Base directory for output
# 3)  Tvan data file path for 30 minute summary in July

user = 'wwieder'
if (user ==  'wwieder') {
  DirBase <- "~/Desktop/Working_files/Niwot/CLM/"
  DirOutBase <- paste0(DirBase,"OBS_SIM_COMP/")
  tvan_data_fp <- "~/Downloads/CLM/datav20200824T1008/data/tvan_forcing_data_precip_mods_both_towers_2007-05-11_2020-08-11.txt"
} else {
  DirBase <- "~/Downloads"
  DirOutBase <- paste0(DirBase,"OBS_SIM_COMP/")
  tvan_data_fp <- "~/Downloads/CLM/datav20200816T1808/data/tvan_forcing_data_precip_mods_both_towers_2007-05-11_2020-08-11.txt"
}


# Simulation Name (for organizing output and naming)
# This is the same as the "case_name" from flow.sim.R
sim_name <- "clm50bgc_NWT_newPHS_lowSLA" #'clm50bgc_NWT_base

#### Input options ####
# Simulation data directory (output from flow.sim.R script)
DirSimIn = paste0(DirBase,'SIM/',sim_name)

# Observation data directory (output from flow.obs.R script)
DirObsIn = paste0(DirBase,'OBS/data')

# What vegetation community are we working with?
vegetation_com <- "SB" # Options: "FF", "DM", "WM", "MM", "SB", NA


##############################################################################
# Static workflow parameters - these are unlikely to change
##############################################################################
DirOut <- paste0(DirOutBase, sim_name)
#Check if directory exists and create if not
if(!dir.exists(DirOut)) dir.create(DirOut, recursive = TRUE)

# simulation file list
sim_file_list <- list.files(DirSimIn, full.names = TRUE)

# observation file list
obs_file_list <- list.files(DirObsIn, full.names = TRUE)


##############################################################################
# Diel plots - fast model/obs timestamp comparison
##############################################################################
# Load CLM simulation data
hlf_hrly_file <- grep("30", sim_file_list)
hlf_hr_flx.clm <- read.table(file = sim_file_list[hlf_hrly_file],
                      sep = "\t", header = TRUE)

# Load Obs data
hlf_hrly_file <- grep("July", obs_file_list)
hlf_hr_flx.obs <- read.table(file = obs_file_list[hlf_hrly_file],
                      sep = "\t", header = TRUE)

# Reformat simulation fluxes
diurnal_flx_vars <- c("RNET", "FSH", "EFLX_LH_TOT", "GPP")

hlf_hr_flx.clm <- hlf_hr_flx.clm %>% 
  filter(veg_com == vegetation_com) %>%
  select(Hour, DoY, year, month, all_of(diurnal_flx_vars)) %>%
  filter(month == 7) %>%
  group_by(Hour) %>%
  summarize_at(all_of(diurnal_flx_vars),
               list(houravg = mean, hoursd = sd), na.rm = TRUE) %>%
  mutate(ObsSim = "Sim")


##############################################################################
# Diel plots
##############################################################################
diel.plot <- hlf_hr_flx.clm %>%
  select(Hour, ends_with("avg")) %>%
  pivot_longer(cols = !Hour, 
               names_to = "Sim_diurnal_flx", 
               values_to = "Sim_value") %>%
  left_join(hlf_hr_flx.obs %>%
              select(Hour, ends_with("avg")) %>%
              pivot_longer(cols = !Hour, 
                           names_to = "Obs_diurnal_flx", 
                           values_to = "Obs_value"), 
            by = c("Hour" = "Hour", "Sim_diurnal_flx" = "Obs_diurnal_flx")) %>%
  rename(diurnal_flx = Sim_diurnal_flx)


diel_plot <- ggplot(data = diel.plot, aes(x = Obs_value, y = Sim_value)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(~diurnal_flx, scales = "free")
diel_plot

##############################################################################
# Load in flux data
##############################################################################
# Load CLM simulation data
diurnal_file <- grep("Diurnal", sim_file_list)
flx.clm <- read.table(file = sim_file_list[diurnal_file],
                      sep = "\t", header = TRUE)

# Load Obs data
diurnal_file <- grep("Diurnal", obs_file_list)
flx.obs <- read.table(file = obs_file_list[diurnal_file],
                      sep = "\t", header = TRUE)


##############################################################################
# Plot flux data
##############################################################################
# Combine flux data for plotting
flx.clm <- flx.clm %>% 
  filter(veg_com == vegetation_com) %>%
  select(all_of(names(flx.obs)))


flx.plot <- bind_rows(flx.clm, flx.obs) %>%
  # reorder months in order of season
  mutate(MonGroup = factor(MonGroup, levels = c("JJA", "MAM", "DJF", "SON")))
  

plot_forcing_var <- function(x) {
  #x <- "RNET"
  plot.df <- flx.plot %>%
    select(MonGroup, Hour, ObsSim, starts_with(x)) %>%
    rename(hourly_mean := !!quo_name(paste0(x, "_houravg")),
           hourly_sd := !!quo_name(paste0(x, "_hoursd")))

  ylabels <- c("GPP" = expression('GPP ('~gC~m^-2~s^-1~')'),
            "FSH" = expression('Sensible Heat Flux ('~W~m^-2~')'),
            "EFLX_LH_TOT" = expression('Latent Heat Flux ('~W~m^-2~')'),
            "RNET" = expression('Net Radiation ('~W~m^-2~')'))

  ylab <- ylabels[x]

  ggplot(plot.df, aes(x = Hour)) +
    geom_ribbon(aes(ymin = hourly_mean - hourly_sd,
                  ymax = hourly_mean + hourly_sd,
                  fill = ObsSim), alpha = 0.5) +
    geom_line(aes(y = hourly_mean, color = ObsSim)) +
    scale_color_manual(values = c("black", "firebrick")) +
    scale_fill_manual(values = c("black", "firebrick")) +
    facet_wrap(~ MonGroup, ncol = 1, strip.position = "left") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(ylab) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "transparent",
                                          color = "transparent"),
          strip.placement = "outsize",
          legend.position = "none")
}

plots <- map(c("RNET", "FSH", "EFLX_LH_TOT","GPP"), ~plot_forcing_var(x = .x))
names(plots) <- c("RNET", "FSH", "EFLX_LH_TOT","GPP")

flx_comp_plot <- cowplot::plot_grid(plotlist = get("plots"), ncol = 4)

cowplot::save_plot(flx_comp_plot, 
                   filename = paste0(DirOut, "/flux_comp_plot_",vegetation_com,".png"),
                   base_height = 6)


##############################################################################
# Load in daily soil moisture, soil temp, and GPP data
##############################################################################
# Load CLM simulation data
daily_file <- grep("Daily", sim_file_list)
daily.clm <- read.table(file = sim_file_list[daily_file],
                      sep = "\t", header = TRUE)
daily.clm <- daily.clm %>% 
  filter(veg_com == vegetation_com) %>%
  select(DoY, ObsSim, veg_com, contains("SOI"), contains("GPP")) 

# Change names to reflect obs names
names(daily.clm) <- sub("TSOI", "soiltemp", names(daily.clm))
names(daily.clm) <- sub("H2OSOI", "soilmoisture", names(daily.clm))
names(daily.clm) <- sub("doyavg", "dailyavg", names(daily.clm))
names(daily.clm) <- sub("doysd", "dailysd", names(daily.clm))

# Load Observational data
daily_file <- grep("Daily", obs_file_list)
daily.obs <- read.table(file = obs_file_list[daily_file],
                        sep = "\t", header = TRUE)

## quick plot of all results
names(daily.obs)
ggplot(daily.obs, aes(x = DoY)) +
  geom_line(aes(y = soilmoisture_upper_avg_dailyavg, color = veg_com))

daily.obs <- daily.obs %>% 
  select(!contains("snow_depth")) %>%
  filter(veg_com == vegetation_com)
names(daily.obs) <- sub("_avg_", "_", names(daily.obs))

##############################################################################
# Plot soil moisture data
##############################################################################

daily.plot <- bind_rows(daily.clm, daily.obs) %>%
  pivot_longer(ends_with("dailyavg"),
               names_to = "MeanMetric",
               values_to = "DailyMean") %>%
  pivot_longer(ends_with("dailysd"),
               names_to = "SDMetric",
               values_to = "DailySD") %>% 
  mutate(MeanMetric = gsub("_dailyavg", "", MeanMetric),
         SDMetric = gsub("_dailysd", "", SDMetric)) %>%
  filter(MeanMetric == SDMetric) %>%
  # change the order of MeanMetric for more intuitive plots
  mutate(MeanMetric = factor(MeanMetric, 
                             levels = c("GPP", "soilmoisture_upper",
                                        "soilmoisture_lower", 
                                        "soiltemp_upper", "soiltemp_lower"))) %>%
  # make a dummy date for easy plotting
  mutate(dummydate = days(DoY) + ymd("2000-01-01"))

soil_moisture_plot <- ggplot(daily.plot, aes(x = dummydate)) +
  geom_ribbon(aes(ymin = DailyMean - DailySD, 
                  ymax = DailyMean + DailySD,
                  fill = ObsSim), alpha = 0.4) +
  geom_line(aes(y = DailyMean, color = ObsSim)) +
  facet_wrap(~MeanMetric, scales = "free_y", ncol = 1) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_color_manual(values = c("black", "firebrick")) +
  scale_fill_manual(values = c("black", "firebrick")) +
  theme_bw() +
  xlab("Day of Year") + ylab("") +
  ggtitle(paste0("Soil properties and GPP for ", vegetation_com, " community"))

ggsave(soil_moisture_plot, 
       file = paste0(DirOut, "/soil_comp_plot_", vegetation_com, ".png"))

##############################################################################
# Load in unsummarized snow depth data
##############################################################################
# load in simulations
clm_file <- grep("Unsummarized", sim_file_list)
all.clm <- read.table(file = sim_file_list[clm_file],
                        sep = "\t", header = TRUE)

# Summarize clm snow depth 
snow_depth.clm <- all.clm %>% 
  select(date, SNOW_DEPTH, veg_com, ObsSim) %>%
  group_by(date, veg_com) %>%
  mutate(avg_snwdp = mean(SNOW_DEPTH, na.rm = TRUE),
         sd_snwdp = sd(SNOW_DEPTH, na.rm = TRUE)) %>%
  select(-SNOW_DEPTH) %>%
  unique()

# Load in observations
snwdp_obs_file <- grep("snow_depth", obs_file_list)

snow_depth.obs <- read.table(file = obs_file_list[snwdp_obs_file],
           sep = "\t", header = TRUE)


# Rename observational data columns to match clm data columns, 
# filter by vegetation community
snow_depth.obs <- snow_depth.obs %>% 
  rename(avg_snwdp = avg_date_depth,
         sd_snwdp = sd_date_depth) %>%
  mutate(ObsSim = "Obs") %>%
  #filter(veg_com == vegetation_com) %>%
  select(-DoY, -data_information, -Year)


# Combine observation and sim data sets
names(snow_depth.obs)
names(snow_depth.clm)


snow_depth.plot <- bind_rows(snow_depth.clm, snow_depth.obs)


snow_depth_plot <- ggplot(snow_depth.plot %>%
         # add a very small number since geom ribbon can't handle widths of 0
         mutate(sd_snwdp = ifelse(sd_snwdp == 0, 0.000000000001, sd_snwdp)),
       aes(x = as.Date(date))) +
  geom_ribbon(aes(ymin = (avg_snwdp - sd_snwdp),
                  ymax = (avg_snwdp + sd_snwdp), 
                  group = ObsSim,
                  fill = ObsSim), alpha = 0.4) +
  geom_line(aes(y = avg_snwdp,
                group = ObsSim,
                color = ObsSim)) +
  facet_wrap(~veg_com, ncol = 1, scales = "free_y") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  scale_color_manual(values = c("black", "firebrick")) +
  scale_fill_manual(values = c("black", "firebrick")) +
  theme_bw() +
  xlab("Day of Year") + ylab("Snow Depth (cm)") +
  ggtitle(paste0("Snow depth"))
  

ggsave(snow_depth_plot, 
       file = paste0(DirOut, "/snow_depth_plot.png"))

print('---- finished plotting ----')
