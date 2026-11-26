rm(list=ls()); graphics.off(); cat("\014")

kerns <- c("gauss","exp")
betas <- c(0,1,3,6)

TABLE <- NULL
for( kern in kerns ) {
  for( b in betas ) {
    
    BOpc <- readRDS(paste0("RESULTS_vanilla_BO_on_pointclouds/results_30_",kern,"_",b,"_11_500_1.RDS"))
    BOpc <- (aggregate(BOpc$y,by=list(BOpc$seed),min))$x
    
    BOfl <- readRDS(paste0("RESULTS_vanilla_BO_on_flows/results_30_",kern,"_",b,"_11_500_1.RDS"))
    BOfl <- (aggregate(BOfl$y,by=list(BOfl$seed),min))$x
    
    PIBO <- readRDS(paste0("RESULTS_PIBO/results_30_",kern,"_",b,"_11_500_1.RDS"))
    PIBO <- (aggregate(PIBO$y,by=list(PIBO$seed),min))$x
    
    TABLE <- rbind(TABLE, data.frame( kern=kern,
                                      beta=b,
                                      test="PIBO vs BO on point clouds",
                                      p.value=round(wilcox.test(PIBO,BOpc,paired=T)$p.value,3) ))
    
    TABLE <- rbind(TABLE, data.frame( kern=kern,
                                      beta=b,
                                      test="PIBO vs BO on flows",
                                      p.value=round(wilcox.test(PIBO,BOfl,paired=T)$p.value,3) ))
    
    TABLE <- rbind(TABLE, data.frame( kern=kern,
                                      beta=b,
                                      test="BO on flows vs BO on point clouds",
                                      p.value=round(wilcox.test(BOfl,BOpc,paired=T)$p.value,3) ))
  }
}

print(TABLE)