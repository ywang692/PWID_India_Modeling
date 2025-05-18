
setup <- function(dat,at) {
  
  if (at == 2){
    dat$`_last_vname` <- tail(sort(data$record_id), 1)
    dat <- set_attr(dat, "active", rep(1,pop))
    dat <- set_attr(dat, "age", data$age)
    dat <- set_attr(dat, "status", status_vec)
    dat <- set_attr(dat, "HIV", get.node.attr(dat$nw[[1]], "HIV"))
      infTime_hiv <- rep(0, pop)
      infTime_hiv[which(hiv_vec == 1)] <- 1
    dat <- set_attr(dat, "infTime_hiv", infTime_hiv)
    dat <- set_attr(dat, "current_inj", rep(1,pop))
    dat <- set_attr(dat, "inj_freq", pmax(data$inj_freq, 1))
    dat <- set_attr(dat, "venues", data$venues)
    dat <- set_epi(dat, "n_current_inj", at, sum(dat$attr$current_inj, na.rm=T))
    
  }
  return (dat)
}


test <- function(dat, at) {

  loc_spec <- function(ids, var){
    out <- rep(NA, length(ids))
    for (i in 1:length(ids)){
      loc <- unlist(dat$nw[[1]][["val"]][[ids[i]]])[venue_list]
      loc <- as.numeric(loc)
      if (sum(loc) > 0){
        out[i] <- max(dat$param$spat_list[[at]][[var]][as.logical(loc)], na.rm=T)
      } else {
        out[i] <- spat_list[[at]][[var]][length(venue_list)+1]
      }
    }
    return(out)
  }
  
  active <- get_attr(dat, "active")
  current_inj <- get_attr(dat, "current_inj")
  alive <- ifelse(is.na(get_attr(dat,"exitTime")), 1, 0)
  
  ## HIV testing
  if (at == 2) {
    diag.status <- rep(0, pop)
    infected <- which(get_attr(dat,"HIV") == 1)
    tested <- infected[as.logical(rbinom(length(infected), 1, prob = loc_spec(infected, "test_ever")))]
    tested <- infected
    diag.status[tested] <- 1
    nTest <- length(tested)
  } else {
    diag.status <- get_attr(dat, "diag.status")
    hiv_status <- get_attr(dat,"HIV")
    infected <- which(hiv_status == 1)

    idsElig <- c(which(alive == 1 & current_inj == 1 & hiv_status == 0), which(alive ==1 & current_inj ==1 & hiv_status == 1 & diag.status == 0))
    nElig <- length(idsElig)

    vecTest <- rbinom(nElig, 1, loc_spec(idsElig, "test_1m"))
    idsTest <- idsElig[as.logical(vecTest)]
    diag.status[idsTest] <- 1
    diag.status[hiv_status==0] <- 0
    nTest <- length(intersect(idsTest, infected))
  }
  dat <- set_attr(dat, "diag.status", diag.status)
  dat <- set_epi(dat, "nTest", at, nTest) 

  return(dat)
}


treatment <- function(dat, at) {

  loc_spec <- function(ids, var){
    out <- rep(NA, length(ids))
    for (i in 1:length(ids)){
      loc <- unlist(dat$nw[[1]][["val"]][[ids[i]]])[venue_list]
      loc <- as.numeric(loc)
      if (sum(loc) > 0){
        out[i] <- max(dat$param$spat_list[[at]][[var]][as.logical(loc)], na.rm=T)
      } else {
        out[i] <- spat_list[[at]][[var]][length(venue_list)+1]
      }
    }
    return(out)
  }
  
  active <- get_attr(dat, "active")
  current_inj <- get_attr(dat, "current_inj")
  test.status <- get_attr(dat, "diag.status")
  hiv.status <- get_attr(dat,"HIV")
  alive <- ifelse(is.na(get_attr(dat,"exitTime")), 1, 0)
  
  ## HIV treatment (ART)
  if (at == 2) {
    suppressed <- data$hiv_suppressed
    suppressed[is.na(suppressed)] <- 0
    dat <- set_attr(dat, "viral.supp", suppressed)
    dat <- set_attr(dat, "art.status", suppressed)
    nArt <- sum(suppressed)
    
  } else {
    art.status <- get_attr(dat, "art.status")
    supp.status <- get_attr(dat, "viral.supp")
    ids.art.on <- which(art.status==1)
    vec_no_art <- rbinom(length(ids.art.on), 1, prob= 1-loc_spec(ids.art.on, "ART_adhere"))
    ids_no_art <- ids.art.on[as.logical(vec_no_art)]
    if (length(ids_no_art) > 0){
      art.status[ids_no_art] <- 0
      idsElig <- which(art.status == 0 & supp.status == 1)
      vec_no_supp <- rbinom(length(idsElig), 1, prob = 1/3)
      if (sum(vec_no_supp) > 0){
        ids_no_supp <- idsElig[as.logical(vec_no_supp)]
        supp.status[ids_no_supp] <- 0
      }
    }
  
    idsElig <- which(alive==1 & hiv.status ==1 & test.status==1 & art.status==0)
    nElig <- length(idsElig)
    if (nElig > 0){
      vec_art <- rbinom(nElig, 1, prob=loc_spec(idsElig, "ART"))
      idsArt <- idsElig[as.logical(vec_art)]
      art.status[idsArt] <- 1
      dat <- set_attr(dat, "art.status", art.status)
    
      vecElig <- alive==1 & hiv.status ==1 & test.status==1 & art.status==1 & supp.status==0
      idsElig <- which(vecElig == T)
      nElig <- length(idsElig)
      if (length(nElig) > 0){
        vec_supp <- rbinom(nElig, 1, prob = 0.92)
        supp.status[idsElig[as.logical(vec_supp)]] <- 1
      }
    }
    dat <- set_attr(dat, "art.status", art.status)
    dat <- set_attr(dat, "viral.supp", supp.status)
  }
  dat <- set_epi(dat, "nArt", at, sum(get_attr(dat, "art.status")))
  dat <- set_epi(dat, "nSupp", at, sum(get_attr(dat, "viral.supp")))

  return(dat)
}


infect <- function(dat, at) {
  
  loc_spec <- function(ids, var){
    out <- rep(NA, length(ids))
    for (i in 1:length(ids)){
      loc <- unlist(dat$nw[[1]][["val"]][[ids[i]]])[venue_list]
      loc <- as.numeric(loc)
      if (sum(loc) > 0){
        out[i] <- max(dat$param$spat_list[[at]][[var]][as.logical(loc)], na.rm=T)
      } else {
        out[i] <- spat_list[[at]][[var]][length(venue_list)+1]
      }
    }
    return(out)
  }
  
  active <- get_attr(dat, "active")
  current_inj <- get_attr(dat, "current_inj")
  status <- get_attr(dat, "status")
  hiv_status <- get_attr(dat, "HIV")
  alive <- is.na(get_attr(dat,"exitTime")) %>% as.numeric()
  idsInf <- which(active == 1 & current_inj == 1 & status == "i" & hiv_status == 1)
  nElig <- length(idsInf)
  
  # To count people who used MOUD/SSP at given timepoint 
  if (at == 2){
   dat$epi$n_moud <- c(0, 0)
   dat$epi$n_ssp <- c(0, 0)
  } else {
   dat$epi$n_moud[at] <- 0
   dat$epi$n_ssp[at] <- 0
  }

  if (nElig > 0) {
    
    #nw <- as.edgelist(dat$nw[[1]]) %>% as.data.frame()
    nw <- data.frame(V1 = unlist(lapply(dat$nw[[1]]$mel, function(x) x$inl)),
                     V2 = unlist(lapply(dat$nw[[1]]$mel, function(x) x$outl)))
    nw$start <- unlist(lapply(dat$nw[[1]]$mel, function(x) as.matrix(x$atl$active[1])))
    nw$end <- unlist(lapply(dat$nw[[1]]$mel, function(x) as.matrix(x$atl$active[2])))
    nw$at <- at >= nw$start & at <= nw$end
    nw <- nw[nw$at, ]
    
    degree <- nrow(nw)*2/sum(active)
    if (at == 2){
      dat$epi$mean_degree <- c(0,degree)
    } else {
      dat$epi$mean_degree[at] <- degree
    }
    
    inj_freq <- get_attr(dat,"inj_freq")[idsInf]
    inj_freq <- pmax(rpois(nElig, inj_freq), 1)
    
    ssp <- rbinom(nElig,1,prob=loc_spec(idsInf, "SSP")) * 0.9
    dat$epi$n_ssp[at] <- sum(ssp > 0)
    moud <- rbinom(nElig,1,prob=loc_spec(idsInf, "MOUD")) * 0.56
    dat$epi$n_moud[at] <- sum(moud > 0)
    
    inj_freq <- inj_freq * (1-moud) * loc_spec(idsInf, "shared_inj")
    
    if (at == 2){
      acute_hiv <- rep(0, nElig)
    } else {
      acute_hiv <- ifelse(get_attr(dat,"infTime_hiv")[idsInf] %in% c(at-2, at-1), 1, 0)
    }
    acute_hiv_increase <- ifelse(acute_hiv == 1, 5.3, 1)
    diagnosed <- get_attr(dat,"diag.status")[idsInf]
    art <- get_attr(dat,"art.status")[idsInf]
    supp <- get_attr(dat,"viral.supp")[idsInf]
    viral_supp_scale <- ifelse(diagnosed + art == 2, 0.728, 0)
    viral_supp_scale[which(diagnosed + art + supp == 3)] <- 0.94
    #hcv_infected <- sample(1:nElig, runif(1, 0.468, 0.709)*nElig, replace=F)
    #viral_supp_scale[hcv_infected] <- viral_supp_scale[hcv_infected] * rbeta(length(hcv_infected), 35, 100)
    prob_syringe_share <- pmax((1 - ssp)* loc_spec(idsInf, "syringe_share"), 0)
    prob_hiv_per_share <- pmax((0.63/100) * acute_hiv_increase * (1 - viral_supp_scale), 0)

    for (ii in 1:nElig){
      index <- idsInf[ii]

      el <- unique(c(nw[which(nw$V1 == index), "V2"], nw[which(nw$V2 == index), "V1"]))
      inj_partners <- intersect(el, which(active == 1))

      if (length(inj_partners)>0){

        ## Injection network
        if (length(inj_partners) > 0){

          n_inj_partners <- length(inj_partners)
          distribution <- get_attr(dat,"inj_freq")[inj_partners]
          
          inj_acts_per_partner <- sample(1:n_inj_partners, size=round(inj_freq[ii]), replace=T, prob=distribution)
          inj_acts_per_partner <- sapply(1:n_inj_partners, function(x) sum(inj_acts_per_partner == x))
          share_syringe_events_per_partner <- sapply(inj_acts_per_partner, function(x) sum(rbinom(x, size = 1, prob = prob_syringe_share[ii])))
          new_hiv <- sapply(1:n_inj_partners, function(x) rbinom(1, size = share_syringe_events_per_partner[x], prob = prob_hiv_per_share[ii]))
          
          idsNewInf_hiv <- inj_partners[which(new_hiv > 0)]

          dat$attr$infTime_hiv[idsNewInf_hiv[as.logical(new_hiv)]] <- at
          dat$attr$HIV[idsNewInf_hiv] <- 1
          dat$attr$status[idsNewInf_hiv] <- "i"
        }
      }
    }
  }
  nInf_hiv <- length(which(dat$attr$HIV ==1 & dat$attr$infTime_hiv == at))
  
  if (at == 2) {
    dat$epi$hiv.si.flow <- c(0, nInf_hiv)
    dat$epi$hiv.i.num <- c(0, sum(active == 1 & dat$attr$HIV ==1))
    dat$epi$s.num <- c(0, sum(active == 1 & dat$attr$status == "s"))
  }
  else {
    dat$epi$hiv.si.flow[at] <- nInf_hiv
    dat$epi$hiv.i.num[at] <- sum(active == 1 & dat$attr$HIV ==1)
    dat$epi$s.num[at] <- sum(active == 1 & dat$attr$status == "s")
  }
  return(dat)
}


aging <- function(dat,at){
  
  alive <- is.na(get_attr(dat,"exitTime"))
  
  dat$attr$age[alive] <- dat$attr$age[alive] + 1/12 
  if (at == 2){
    dat$epi$meanAge <- c(NA_real_, mean(dat$attr$age[alive], na.rm=T))
  } else {
    dat$epi$meanAge[at] <- mean(dat$attr$age[alive], na.rm=T)
  }
  return(dat)
}


dfunc <- function(dat,at){
  
  alive <- get_attr(dat,"exitTime") %>% sapply(function(x) ifelse(is.na(x), 1, 0))
  idsElig <- which(alive ==1)
  nElig <- length(idsElig)
  nDeaths <- 0
  if (nElig > 0) {
    hiv <- as.logical(dat$attr$HIV[idsElig] == 1)
    
    death.rates <- runif(nElig, min=0.0069059, max=0.0116738)
    # death.rates <- death.rates / (0.37 * 1.17 + (1 - 0.37)) #first set everyone to have HIV- mortality
    # death.rates[hiv] <- death.rates[hiv] * 1.17 #then multiply mortality of HIV+ indiv
    vecDeaths <- which(rbinom(nElig, 1, death.rates) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)
    if (nDeaths > 0){
      dat$attr$active[idsDeaths] <- 0
      dat$attr$exitTime[idsDeaths] <- at
      deactivate.vertices(dat$nw[[1]], v = idsDeaths, onset = at, terminus = Inf, deactivate.edges = TRUE)
    }
  }
  
  if (at == 2) {
    dat$epi$d.flow <- c(0, nDeaths)
  } else {
    dat$epi$d.flow[at] <- nDeaths
  }
  return(dat)
}


bfunc <- function(dat, at){

  n_nw <- network.size(dat$nw[[1]])
  
  alive <- get_attr(dat,"exitTime") %>% sapply(function(x) ifelse(is.na(x), 1, 0))
  active <- get_attr(dat, "active")
  
  numNeeded <- pop - sum(dat$attr$active == 1)
  
  if (numNeeded > 0){
    nBirths <- rpois(1, numNeeded)
  } else {
    nBirths <- 0
  }
  if (nBirths > 0){
    last_pid <- dat$`_last_unique_id`
    last_vname <- dat$`_last_vname`
    last_vname_num <- as.numeric(substring(last_vname, 3,6))
    new_vname_num <- seq(1,nBirths) + last_vname_num
    new_vname <- paste0("DC", new_vname_num, sep="")
    new_pid <- seq(1,nBirths) + last_pid
    dat[["_last_unique_id"]] <- tail(new_pid, 1)
    dat[["_last_vname"]] <- tail(new_vname, 1)
    dat$nw[[1]] <- add.vertices.networkDynamic(dat$nw[[1]], nv=nBirths, vertex.pid = as.character(new_pid))
    network.vertex.names(dat$nw[[1]])[new_pid] <- new_vname
    newNodes <- (n_nw + 1) : (n_nw + nBirths)
    dat$attr$unique_id <- c(dat$attr$unique_id, new_pid)
    
    dat$attr$active <- c(dat$attr$active, rep(1,nBirths))
      hiv_status <- rbinom(nBirths, 1, prob=0.093) #naco.gov.in(2.1%)
    dat$attr$HIV <- dat$attr$HIV <- c(dat$attr$HIV, hiv_status)
    dat$attr$infTime_hiv <- c(dat$attr$infTime_hiv, rep(0, nBirths))
      status <- rep("s", nBirths)
      status[which(hiv_status == 1)] <- "i"
    dat$attr$status <- c(dat$attr$status, status)
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
    dat$attr$current_inj <- c(dat$attr$current_inj, rep(1, nBirths))
      unborn_age <- pmax(rgamma(nBirths, shape=age_1st_inj_dist[1], rate=age_1st_inj_dist[2]), 18)
    dat$attr$age <- c(dat$attr$age, unborn_age)
      new_freq <- pmax(rgamma(nBirths, shape=inj30d_dist[1], rate=inj30d_dist[2]), 1) %>% round()
    dat$attr$inj_freq <- c(dat$attr$inj_freq, new_freq)
      new_venues <- rgamma(nBirths, shape=venues_dist[1], rate=venues_dist[2]) %>% round() %>% pmin(.,11)
    dat$attr$venues <- c(dat$attr$venues, new_venues)
    
    dat$attr$diag.status <- c(dat$attr$diag.status, rep(0, nBirths))
    dat$attr$art.status <- c(dat$attr$art.status, rep(0, nBirths))
    dat$attr$viral.supp <- c(dat$attr$viral.supp, rep(0, nBirths))
    
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nBirths))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nBirths))
    
    dat$nw[[1]] <- activate.vertices(dat$nw[[1]], v = newNodes, onset = at, terminus = Inf)
  }
  
  if (at == 2){
    dat$epi$b.flow <- c(0,nBirths)
  } else {
    dat$epi$b.flow[at] <- nBirths
  }
  
  nActive <- length(which(dat$attr$active[at] == 1))
  if (at == 2){
    dat$epi$nActive <- c(0,nActive)
  } else {
    dat$epi$nActive[at] <- nActive
  }
  
  current <- which(dat$attr$active == 1 & dat$attr$current_inj == 1)
  dat$epi$n_current_inj[at] <- length(current)
  
  return(dat)
}


nwup <- function(dat, at){
  
  pop_size <- length(dat$attr$active)
  pop_size_pre <- length(dat$attr[["spat_jb1___40"]])
  new <- (pop_size_pre+1) : pop_size
  
  ## Attributes
  status <- get_attr(dat, "status")
  active <- get_attr(dat, "active")
  entrTime <- get_attr(dat, "entrTime")
  exitTime <- get_attr(dat, "exitTime")
  
  ## Controls
  tergmLite <- get_control(dat, "tergmLite")
  resimulate.network <- get_control(dat, "resimulate.network")
  
  ## Vital Dynamics
  arrivals <- which(entrTime == at)
  departures <- which(exitTime == at)
  
  nArrivals <- length(arrivals)
  if (nArrivals > 0) {
    
    nwterms <- dat$temp$nwterms
    if (!is.null(nwterms)) {
      #curr.tab <- get_attr_prop(dat, nwterms)
      #dat <- auto_update_attr(dat, arrivals, curr.tab)
      
      for (i in 1:length(venues)){
        vname <- venues[i]
        dat$attr[[vname]][new] <- 0
      }
      for (i in new){
        v_count <- dat$attr$venues[i]
        if (v_count > 0){
          new_v <- sample(spat_var$loc[1:length(venues)], v_count, prob = spat_var$pop[1:length(venues)], replace=F)
          for (j in 1:length(new_v)){
            dat$attr[[new_v[j]]][i] <- 1
          }
        }
      }
    }
    if (length(unique(sapply(dat$attr, length))) != 1) {
      stop("Attribute list of unequal length. Check arrivals.net module.\n",
           print(cbind(sapply(get_attr_list(dat), length))))
    }
  }
  
  ## Copy static attributes to network object
  if (tergmLite == FALSE && resimulate.network == TRUE) {
    dat <- copy_datattr_to_nwattr(dat)
  }
  
  ## Update temporally extended disease status
  if (tergmLite == FALSE) {
    for (network in seq_len(dat$num.nw)) {
      dat$nw[[network]] <- activate.vertex.attribute(dat$nw[[network]],
                                                     prefix = "testatus",
                                                     value = status,
                                                     onset = at,
                                                     terminus = Inf)
    }
  }
  return(dat)
}


vb <- function (x, type, s = 1, at = 2) {
  
  if (type == "startup") {
    if (x$verbose == TRUE) {
      cat("\nStarting Network Simulation...")
    }
  }
  if (type == "progress") {
    if (x$control$verbose == TRUE) {
      if (x$control$verbose.int == 0 && at == x$control$nsteps) {
        cat("\nSim = ", s, "/", x$control$nsims, sep = "")
      }
      if (x$control$verbose.int > 0 && (at%%x$control$verbose.int == 
                                        0)) {
        cat("\f")
        cat("\nEpidemic Simulation")
        cat("\n----------------------------")
        cat("\nSimulation: ", s, "/", x$control$nsims, 
            sep = "")
        cat("\nTimestep: ", at, "/", x$control$nsteps, 
            sep = "")
        active <- x$attr$active
        cat("\nPopulation Size:", sum(active == 1))
        cat("\nHIV Prevalence:", round(x$epi$hiv.i.num[at]/sum(active == 1), 3))
        cat("\n----------------------------")
      }
    }
  }
}


