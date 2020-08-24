
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

#Install packages from github repos
# devtools::install_github(c("NEONScience/eddy4R/pack/eddy4R.base", "NEONScience/NEON-utilities/neonUtilities"))
# 
#Setup Environment
options(stringsAsFactors = F)

##############################################################################
#Workflow parameters
##############################################################################

#### Output Options ####

# Base directory for output
DirOutBase <- paste0("~/Downloads/OBS_SIM_COMP/")

# Simulation Name (for organizing output and naming)
sim_name <- "2000datm_CLM50bgc_nwt_DM.clm2.h0.2008-01-01-00000"

#### Input options ####
# Simulation data directory (output from flow.sim.R script)
DirSimIn = "~/Downloads/SIM/data/2000datm_CLM50bgc_nwt_DM.clm2.h0.2008-01-01-00000"

# Observation data directory (output from flow.obs.R script)
DirObsIn = "~/Downloads/OBS/data"

# What vegetation community are we working with
vegetation_com <- "FF"
##############################################################################
# Static workflow parameters - these are unlikely to change
##############################################################################
DirOut <- paste0(DirOutBase, sim_name)
#Check if directory exists and create if not
if(!dir.exists(DirOutBase)) dir.create(DirOut, recursive = TRUE)

# simulation file list
sim_file_list <- list.files(DirSimIn, full.names = TRUE)

# observation file list
obs_file_list <- list.files(DirObsIn, full.names = TRUE)
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
  select(all_of(names(flx.obs)))

flx.plot <- bind_rows(flx.clm, flx.obs)

plot_forcing_var <- function(x) {
  #x <- "RNET"
  plot.df <- flx.plot %>%
    select(MonGroup, Hour, ObsSim, starts_with(x)) %>%
    rename(hourly_mean := !!quo_name(paste0(x, "_houravg")),
           hourly_sd := !!quo_name(paste0(x, "_hoursd")))

  ylabels <- c("GPP" = expression('GPP ('~gC~m^-2~s^-2~')'),
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
                   filename = paste0(DirOut, "/flux_comp_plot.png"),
                   base_height = 6)


##############################################################################
# Load in daily soil moisture, soil temp, and GPP data
##############################################################################
# Load CLM simulation data
daily_file <- grep("Daily", sim_file_list)
daily.clm <- read.table(file = sim_file_list[daily_file],
                      sep = "\t", header = TRUE)
daily.clm <- daily.clm %>% 
  select(DoY, ObsSim, veg_com, contains("SOI"), contains("GPP")) 

# Change names to reflect obs names
names(daily.clm) <- sub("TSOI", "soiltemp", names(daily.clm))
names(daily.clm) <- sub("H2OSOI", "soilmoisture", names(daily.clm))
names(daily.clm) <- sub("_5_", "_lower_", names(daily.clm))
names(daily.clm) <- sub("_3_|_2_", "_upper_", names(daily.clm))
names(daily.clm) <- sub("doyavg", "dailyavg", names(daily.clm))
names(daily.clm) <- sub("doysd", "dailysd", names(daily.clm))

# Load Observational data
daily_file <- grep("Daily", obs_file_list)
daily.obs <- read.table(file = obs_file_list[daily_file],
                        sep = "\t", header = TRUE)

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
       file = paste0(DirOut, "/soil_comp_plot.png"))


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
  group_by(date) %>%
  mutate(avg_snwdp = mean(SNOW_DEPTH, na.rm = TRUE),
         sd_snwdp = sd(SNOW_DEPTH, na.rm = TRUE))

# Load in observations
snwdp_obs_file <- grep("snow_depth", obs_file_list)

snow_depth.obs <- read.table(file = obs_file_list[snwdp_obs_file],
           sep = "\t", header = TRUE)

snow_depth.obs %>% names()
  rename()


snw_dpth.clm <- snw_dpth.clm %>% 
  select(DoY, ObsSim, veg_com, contains("SOI"), contains("GPP")) 


