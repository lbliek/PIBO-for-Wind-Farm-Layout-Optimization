rm(list=ls()); graphics.off(); cat("\014")

folders <- c("RESULTS_PIBO", "RESULTS_vanilla_BO_on_flows","RESULTS_vanilla_BO_on_pointclouds")

for( folder in folders ) {
  
  if( folder=="RESULTS_PIBO" ) {
    cat(" [ PIDO ]\n")
  } else {
    if( folder=="RESULTS_vanilla_BO_on_flows" ) {
      cat(" [vannila BO on Flows]\n")
    } else {
      cat(" [vanilla BO on point-clouds (Physical Space) ]\n")
    }
  }
  
  allFiles <- list.files(folder)
  
  TABLE <- NULL
  for( f in allFiles ) {
    res <- readRDS( paste0(folder,"/",f) )
    
    aggr <- aggregate(res$y,by=list(res$seed),min)
    
    info <- unlist(strsplit(f,"_",fixed=F))
    nSeeds <- as.numeric(info[2])
    kernel <- info[3]
    lcb.beta <- as.numeric(info[4])
    n0 <- as.numeric(info[5])
    N <- as.numeric(info[6])
    
    TABLE <- rbind( TABLE, data.frame( nSeeds=nSeeds,
                                       kernel=kernel,
                                       lcb.beta=lcb.beta,
                                       n0=n0,
                                       N=N,
                                       bs_avg=round(mean(aggr$x),2),
                                       bs_sd=round(sd(aggr$x),2),
                                       bs_med=round(median(aggr$x),2),
                                       bs_min=round(min(aggr$x),2),
                                       bs_max=round(max(aggr$x),2),
                                       stringsAsFactors=F ))
  }
  
  print(TABLE)
  cat("\n\n")
}


