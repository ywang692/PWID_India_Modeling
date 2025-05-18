library(tidyverse)
library(fitdistrplus)

# Load data ---------------------------------------------------------------

data <- readRDS("data/cleaned/data_indiv.rds")
el <- readRDS("data/raw/edgelist_synthetic.rds")
spat_var <- readRDS("data/cleaned/data_venue.rds")
spat_var <- spat_var[order(spat_var$pop, decreasing=T), ]

venues <- spat_var$loc[1:(nrow(spat_var)-1)]
venue_list <- venues

el <- el[order(el$Source),]
data <- data[c(data$record_id) %in% unlist(el),]


# Params ------------------------------------------------------------------

ids <- unique(c(el$Source, el$Target))
pop <- n <-  nrow(data)
n_years <- 5 #run simulations for 5 years including burn-in period
time_steps <- n_years * 12 #every time step is 1 month

# Create distributions for new arrivals into model simulations
age_1st_inj_dist <- fitdistr(data$age_1st_inj[!is.na(data$age_1st_inj)], "gamma") %>% coefficients()
inj30d_dist <- fitdistr(pmax(data$inj_freq[!is.na(data$inj_freq)],1), "gamma") %>% coefficients()
venues_dist <- fitdistr(pmin(pmax(data$venues, 0.1),8), "gamma") %>% coefficients()

