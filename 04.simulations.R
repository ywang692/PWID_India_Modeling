source("01.param.R")
source("02.net_gen.R")
source("03.API_extensions.R")

hiv_vec <- get.node.attr(est1$newnetwork, "HIV") %>% as.numeric()

status_vec <- rep("s", pop)
status_vec[which(hiv_vec==1)] <- "i"

init <- init.net(status.vector = status_vec)

module.order <- c("setup.FUN",
                  "resim_nets.FUN",
                  "test.FUN",
                  "treatment.FUN",
                  "infection.FUN",
                  "aging.FUN",
                  "departures.FUN",
                  "arrivals.FUN",
                  "nwupdate.FUN",
                  "prevalence.FUN")

control <- control.net(type=NULL, nsims=1, nsteps=time_steps,
                       ncores=1,
                       setup.FUN = setup,
                       test.FUN = test,
                       treatment.FUN = treatment,
                       infection.FUN = infect,
                       aging.FUN = aging,
                       departures.FUN = dfunc,
                       arrivals.FUN = bfunc,
                       nwupdate.FUN = nwup,
                       verbose.FUN = vb,
                       resimulate.network=T,
                       module.order = module.order,
                       verbose = T,
                       save.trans=F)


ssp_cov <- c(0.55, 0.41, 0.28)
moud_cov <- c(0.50, 0.38, 0.25)
test_cov <- c(0.95, 0.71, 0.48)
art_cov <- c(0.95, 0.71, 0.48)

coverage <- c(1, 2, 3)
cov_df <- expand.grid(coverage, coverage, coverage, coverage)
cov_df <- cov_df[apply(cov_df, 1,
                       function(x) x[1] <= x[2] && x[2] <= x[3] && x[3] <= x[4]), ]
colnames(cov_df) <- c("SSP", "MOUD", "test_1m", "ART")

cov_df$SSP <- sapply(cov_df$SSP, function(x) ssp_cov[as.numeric(x)])
cov_df$MOUD <- sapply(cov_df$MOUD, function(x) moud_cov[as.numeric(x)])
cov_df$test_1m <- sapply(cov_df$test_1m, function(x) test_cov[as.numeric(x)])
cov_df$ART <- sapply(cov_df$ART, function(x) art_cov[as.numeric(x)])


coef.diss.list <- list()
coef.diss.list[1:time_steps] <- list(coef.diss)

by_pop <- c("spat_jb1___40", "spat_jb1___39", "spat_jb1___38", "spat_jb1___37", "spat_jp5___105",
            "spat_jb1___41", "spat_jb2___44", "spat_jb2___48", "spat_jp3___165", "spat_nkp___52",  
            "spat_jp5___104", "spat_jb2___47")
by_cluster <- c("spat_jb1___40", "spat_jp5___105", "spat_dk2___67", "spat_jb1___39", "spat_jp3___165", 
                "spat_dk1___65", "spat_jb1___38", "spat_jp5___104", "spat_dk2___68", "spat_jb1___37",  
                "spat_jp5___76", "spat_dk1___66")

v <- match(by_cluster, spat_var$loc)

cov <- 1
vs <- 4

spat_var_new <- spat_var
spat_var_new[v[1:vs], "SSP"] <- cov_df[cov, "SSP"]
spat_var_new[v[1:vs], "MOUD"] <- cov_df[cov, "MOUD"]
spat_var_new[v[1:vs], "test_1m"] <- cov_df[cov, "test_1m"]
spat_var_new[v[1:vs], "ART"] <- cov_df[cov, "ART"]

spat_list <- list()
spat_list[1:40] <- list(spat_var)
spat_list[41:time_steps] <- list(spat_var_new)

param <- param.net(spat_list = spat_list)

sim <- netsim(est1, param, init, control)


