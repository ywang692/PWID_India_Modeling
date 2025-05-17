################################################################################
# This R file contains codes to extract individual-level and injection venue-level parameters for HIV transmission model. 
# The synthetic datasets are in the forms of .csv and .rds. 
# Columns in the synthetic dataset resemble the actual raw dataset and need to be cleaned and compiled. 
# The code below generates two datasets in the folder ~/data/cleaned:
# data_indiv = characteristics of hypothetical participants
# data_venue = characteristics of injection venues by compiling all participants who frequented an injection venue
# 
# Variable names correspond to these exact questions from the survey:
# record_id: unique participant ID
# redcap_event_name: specifies baseline dataset
# age: age calculated from participants' date of birth by the survey date OR directly from the question 'How old are you?'
# age_1st_inj: 'How old were you when you first injected drugs for non-medicinal purposes?'
# hiv: HIV-1 antibody results (1 = reactive, 0 = non-reactive)
# hiv_suppressed: HIV-1 detectable viral load (1 = '≥ 150 copies/mL', 0 = '< 150 copies/mL')
# hiv_aware: 'What were the results of your last HIV test?'
# inj_days_6mo: 'How many days total in the last 6 months did you inject any drugs?'
# inj_times_per_day: 'On days you injected, how many times a day did you inject?'
# MOUD_6mo: 'Have you used an opiate substitution program in the prior 6 months?'
# MOUD_days: 'How many days did you use an opioid substitution program in the prior 6 months?'
# SSP_6mo: 'Have you used a needle syringe exchange program in the prior 6 months?'
# SSP_days: 'How many days did you use the needle exchange program in the last 6 months?'
# ART_ever: 'Have you ever taken antiretroviral medications (tablets to control your HIV infection)?'
# ART_adhere: 'In the past 30 days, did you take any HIV medication?'
# test_ever: 'Have you ever been tested for HIV?'
# test_1mo: 'What is the best estimate of when you were last tested for HIV?' (1 = 'within the last month', 0 = more than 1 month)
# net_X_share (X from 1 to 20): 'How many times have you injected with XX in the prior 30 days? By injected I mean that you were in the same place when you injected and you may or may not have shared injection equipment with them.'
# net_X_syringe_days (X from 1 to 20): 'How many days have you shared needles with XX out of the past 30 days?'
# net_X_syringe_times (X from 1 to 20): 'How many times a day did you share needles with XX in the last 30 days?'
# spat_xxx_xx (identify from a map): 'Take a moment to think about the places where you may have injected in the past 6 months. It could be your house, your dealer s house, at a friend s place, 
# or a public space like a street corner or park or public toilet. If it helps, think back to the last day you injected and where your day started. Over the past 6 months, in which neighborhoods have you injected?'

dat <- readRDS("data/raw/data_synthetic.rds")


# Individual-level dataset ------------------------------------------------

data_indiv <- data.frame(record_id = dat$record_id)
data_indiv$age <- round(dat$age)
data_indiv$age_1st_inj <- round(dat$age_1st_inj)
data_indiv$hiv <- dat$hiv
data_indiv$hiv_suppressed <- ifelse(data_indiv$hiv == 0, NA, as.numeric(dat$hivvl == 0))
data_indiv$hiv_aware <- dat$hiv_aware
data_indiv$inj_freq <- round(dat$inj_days_6mo*dat$inj_times_per_day / 6) #calculate monthly injection frequency
data_indiv$MOUD <- ifelse(dat$MOUD_6mo == 0, 0, dat$MOUD_days)
data_indiv$SSP <- ifelse(dat$SSP_6mo == 0, 0, dat$SSP_days)
data_indiv$ART <- dat$ART_ever
data_indiv$ART_adhere <- dat$ART_adhere
data_indiv$test_ever <- dat$test_ever
data_indiv$test_1m <- dat$test_1m
shared_inj <- dat[ , paste0("net_", 1:20, "_share")] 
shared_inj <- rowSums(shared_inj, na.rm=T) #compile to ego-level total injections with others
data_indiv$shared_inj <- pmin(shared_inj/data_indiv$inj_freq, 1) #calculate proportion of injections done with others
syringe_days <- dat[ , paste0("net_", 1:20, "_syringe_days")] 
syringe_days <- rowSums(syringe_days, na.rm=T)
syringe_times <- dat[ , paste0("net_", 1:20, "_syringe_times")]
syringe_times <- rowSums(syringe_times, na.rm=T)
data_indiv$syringe_share <- pmin(syringe_days * syringe_times / shared_inj, 1) #calculate proportion of injections with others that were done with shared syringes
data_indiv$venues <- rowSums(dat[,grepl("spat_", colnames(dat))], na.rm=T)
data_indiv <- cbind(data_indiv, dat[,grepl("spat_", colnames(dat))])


# Venue-level dataset -----------------------------------------------------

locs <- data_indiv[ ,grepl("spat_", colnames(data_indiv))]
data_venue <- data.frame(loc = colnames(locs),
                     pop = colSums(locs))
data_venue <- rbind(data_venue, c("other", NA))

for (ii in 1:nrow(data_venue)){
  
  loc <- data_venue$loc[ii]
  
  if (loc == "other"){
    tmp <- data_indiv[which(data_indiv$venues == 0), ]
  } else {
    tmp <- data_indiv[which(data_indiv[,loc] == 1), ]
  }
  
  data_venue$locs_per_indiv[ii] <- mean(tmp$venues)
  data_venue$age[ii] <- mean(tmp$age, na.rm=T) 
  data_venue$age_1st_inj[ii] <- mean(tmp$age_1st_inj, na.rm=T)
  data_venue$HIV_active[ii] <- length(which(tmp$hiv ==1)) / nrow(tmp)
  data_venue$HIV_aware[ii] <- length(which(tmp$hiv_aware ==1)) / length(which(tmp$hiv ==1))
  data_venue$HIV_suppressed[ii] <- length(which(tmp$hiv_suppressed == 1)) / length(which(tmp$hiv_aware == 1))
  data_venue$inj_freq[ii] <- mean(tmp$inj_freq, na.rm=T) 
  data_venue$shared_inj[ii] <- mean(tmp$shared_inj, na.rm=T)
  data_venue$syringe_share[ii] <- mean(tmp$syringe_share, na.rm=T) 
  data_venue$MOUD[ii] <- mean(tmp$MOUD/180, na.rm=T) 
  data_venue$SSP[ii] <- mean(tmp$SSP/180, na.rm=T) 
  data_venue$test_ever[ii] <- length(which(tmp$test_ever ==T)) / nrow(tmp)
  data_venue$test_1m[ii] <- length(which(tmp$test_1m == 1)) / length(which(tmp$test_ever == 1))
  data_venue$ART[ii] <- length(which(tmp$ART == 1)) / length(which(tmp$hiv_aware == 1))
  data_venue$ART_adhere[ii] <- length(which(tmp$ART_adhere == 1)) / length(which(tmp$ART == 1))
}

for (ii in 4:ncol(data_venue)){
  col <- data_venue[ ,ii]
  col[is.nan(col)] <- 0
  data_venue[ ,ii] <- col
}


