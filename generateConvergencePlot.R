rm(list=ls()); graphics.off(); cat("\014")

library(mvtnorm)
m <- 5        # number of turbines
rho <- 0.1512 # from 2 * rotor diam / farm length = (2*126)/(333.33*5)

PIBO_res <- readRDS("RESULTS_PIBO/results_30_exp_6_11_500_1.RDS")
BOfl_res <- readRDS("RESULTS_vanilla_BO_on_flows/results_30_exp_6_11_500_1.RDS")
BOpc_res <- readRDS("RESULTS_vanilla_BO_on_pointclouds/results_30_exp_6_11_500_1.RDS")
hyperopt <- read.table("TPE_RS_WWsorted_surrogate_iters.csv",
                     header=T, sep="," )
RS_res <- hyperopt[hyperopt$approach=="hyperopt/randomsearch",]
TPE_res <- hyperopt[hyperopt$approach=="hyperopt/tpe",]



#****************************************************************************
# Convergence plots
#****************************************************************************

# revert to maximization problem...
PIBO_res$y <- -PIBO_res$y
BOfl_res$y <- -BOfl_res$y
BOpc_res$y <- -BOpc_res$y
RS_res$iter_best_fitness <- -RS_res$iter_best_fitness
TPE_res$iter_best_fitness <- -TPE_res$iter_best_fitness


PIBObs <- BOflbs <- BOpcbs <- NULL
for( seed in sort(unique(PIBO_res$seed)) ) {
  
  ix1 <- which(PIBO_res$seed==seed & PIBO_res$iter==0)
  ix2 <- PIBO_res$seed==seed & PIBO_res$iter!=0
  PIBObs <- rbind( PIBObs,
                   cummax( c( max(PIBO_res$y[ix1]), PIBO_res$y[ix2] ) ) )
  
  ix1 <- which(BOfl_res$seed==seed & BOfl_res$iter==0)
  ix2 <- BOfl_res$seed==seed & BOfl_res$iter!=0
  BOflbs <- rbind( BOflbs,
                   cummax( c( max(BOfl_res$y[BOfl_res$seed==seed & BOfl_res$iter==0]),
                              BOfl_res$y[BOfl_res$seed==seed & BOfl_res$iter!=0] ) ) )
  
  ix1 <- which(BOpc_res$seed==seed & BOpc_res$iter==0)
  ix2 <- BOpc_res$seed==seed & BOpc_res$iter!=0
  BOpcbs <- rbind( BOpcbs,
                   cummax( c( max(BOpc_res$y[BOpc_res$seed==seed & BOpc_res$iter==0]),
                              BOpc_res$y[BOpc_res$seed==seed & BOpc_res$iter!=0] ) ) )
}

PIBOm <- apply(PIBObs,2,mean); PIBOlo <- PIBOm-apply(PIBObs,2,sd); PIBOup <- PIBOm+apply(PIBObs,2,sd)
BOflm <- apply(BOflbs,2,mean); BOfllo <- BOflm-apply(BOflbs,2,sd); BOflup <- BOflm+apply(BOflbs,2,sd)
BOpcm <- apply(BOpcbs,2,mean); BOpclo <- BOpcm-apply(BOpcbs,2,sd); BOpcup <- BOpcm+apply(BOpcbs,2,sd)


N <- nrow(PIBO_res)/length(unique(PIBO_res$seed))

RSbs <- TPEbs <- NULL
for( i in 0:4 ) {
  RSbs <- rbind( RSbs, RS_res$iter_best_fitness[(i*N+1):((i+1)*N)][-c(1:(2*m+1))] )
  TPEbs <- rbind( TPEbs, TPE_res$iter_best_fitness[(i*N+1):((i+1)*N)][-c(1:(2*m+1))] )
}
RSm <- apply(RSbs,2,mean); RSlo <- RSm - apply(RSbs,2,sd); RSup <- RSm + apply(RSbs,2,sd)
TPEm <- apply(TPEbs,2,mean); TPElo <- TPEm - apply(TPEbs,2,sd); TPEup <- TPEm + apply(TPEbs,2,sd)





par(mar=c(4.1,5.1,1.1,1.1))
plot( NA, NA, xlim=c(0,length(PIBOm)-1),
      ylim=c(54,71),
      xlab="iterations", ylab="Generated Power", cex.lab=1.5, cex.axis=1.5 )
polygon( c(1:length(TPEm)-1,rev(1:length(TPEm)-1)), c(TPElo,rev(TPEup)), col=adjustcolor("pink",alpha.f=0.4), border="pink" )         
polygon( c(1:length(BOpcm)-1,rev(1:length(BOpcm)-1)), c(BOpclo,rev(BOpcup)), col=adjustcolor("orange",alpha.f=0.4), border="orange" )
polygon( c(1:length(BOflm)-1,rev(1:length(BOflm)-1)), c(BOfllo,rev(BOflup)), col=adjustcolor("green4",alpha.f=0.2), border="green4" )
polygon( c(1:length(PIBOm)-1,rev(1:length(PIBOm)-1)), c(PIBOlo,rev(PIBOup)), col=adjustcolor("deepskyblue",alpha.f=0.4), border="deepskyblue" )

lines( 1:length(TPEm)-1, TPEm, lwd=3, col="firebrick" )
lines( 1:length(BOpcm)-1, BOpcm, lwd=3, col="orange4" )
lines( 1:length(PIBOm)-1, BOflm, lwd=3, col="green4" )
lines( 1:length(PIBOm)-1, PIBOm, lwd=3, col="blue" )

ixs <- seq(1,length(RSm),by=20)
points( ixs-1, BOpcm[ixs], pch=22, cex=1.5, bg="orange" )
points( ixs-1, TPEm[ixs], pch=25, cex=1.2, bg="red" )
points( ixs-1, BOflm[ixs], pch=23, cex=1.4, bg="green4" )
points( ixs-1, PIBOm[ixs], pch=21, cex=1.4, bg="blue" )


points( length(PIBOm)-1, RSm[length(RSm)], pch=8, cex=1.5, lwd=3 )
lines( rep(length(PIBOm)-1,2), c(RSlo[length(PIBOm)-1],RSup[length(PIBOm)-1]),lwd=2)  
lines( length(PIBOm)-1+c(-5,5), rep(RSlo[length(PIBOm)-1],2),lwd=2)  
lines( length(PIBOm)-1+c(-5,5), rep(RSup[length(PIBOm)-1],2),lwd=2)  


legend("bottomright",
       legend=c("PIBO",
                "vanilla BO on flows",
                "vanilla BO on point-clouds",
                "TPE","Random Search"),
       pt.lwd=c(rep(1,4),2),
       pt.bg=c("blue","green4","orange","red","black"),
       pch=c(21,23,22,25,8), cex=1.5, 
       pt.cex=c(1.5,1.5,1.5,1.3,1.5) )






#****************************************************************************


PIBO_tt <- aggregate(PIBO_res$trnTime,by=list(PIBO_res$seed),sum,na.rm=T)
PIBO_at <- aggregate(PIBO_res$acqTime,by=list(PIBO_res$seed),sum,na.rm=T)
PIBO_time <- PIBO_tt$x + PIBO_at$x


BOfl_tt <- aggregate(BOfl_res$trnTime,by=list(BOfl_res$seed),sum,na.rm=T)
BOfl_at <- aggregate(BOfl_res$acqTime,by=list(BOfl_res$seed),sum,na.rm=T)
BOfl_time <- BOfl_tt$x + BOfl_at$x


BOpc_tt <- aggregate(BOpc_res$trnTime,by=list(BOpc_res$seed),sum,na.rm=T)
BOpc_at <- aggregate(BOpc_res$acqTime,by=list(BOpc_res$seed),sum,na.rm=T)
BOpc_time <- BOpc_tt$x + BOpc_at$x


cat("> PIBO - time [secs]:",round(mean(PIBO_time),3),"+/-",round(sd(PIBO_time),3),"\n")
cat("> vanilla BO on flows - time [secs]:",round(mean(BOfl_time),3),"+/-",round(sd(BOfl_time),3),"\n")
cat("> vanilla BO on point-clouds - time [secs]:",round(mean(BOpc_time),3),"+/-",round(sd(BOpc_time),3),"\n")
