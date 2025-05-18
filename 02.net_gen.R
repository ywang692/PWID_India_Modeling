library(EpiModel)

g <- as.network(el, directed=FALSE) #generate network from synthetic edgelist
data <- data[match(get.node.attr(g, "vertex.names"), data$record_id), ] #make sure node indices match participant IDs

g %v% "venues" <- c(data$venues)
g %v% "HIV" <- c(data$hiv)

for (ii in 1:length(venues)){
  g %v% venues[ii] <- c(data[ ,venues[ii]])
}

popular_venues <- venues[1:20]

terms <- paste0("nodematch('", popular_venues, "')")
terms <- c("edges", "nodecov('venues')", "concurrent", "degree(0)", 
           "nodefactor('HIV')", "nodematch('HIV')", terms)
formation <- as.formula(paste("~", paste(terms, collapse = " + ")))
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 18)


est1 <- netest(g, formation, target.stats=NULL, coef.diss)