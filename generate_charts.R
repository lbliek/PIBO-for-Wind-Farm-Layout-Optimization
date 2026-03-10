rm(list=ls()); graphics.off(); cat("\014")

library(mvtnorm)

m <- 5        # number of turbines
rho <- 0.1512 # from 2 * rotor diam / farm length = (2*126)/(333.33*5)

ref_mean <- c(-1,-1)      # bi-variate mean of the reference point-cloud
ref_covm <- 0.01*diag(2)  # covariance matrix of the reference point-cloud  

# Generating the reference point-cloud
refOk <- F; set.seed(42) # this seed is used only for the generation of P_, which will be the same over all the independent runs!
while( !refOk ) {
  P_ <- rmvnorm( m, ref_mean, ref_covm )
  aux <- as.matrix(dist(P_))
  diag(aux) <- max(aux)
  refOk <- (min(aux)<=rho)
}



kern <- "exp"; xi=6 #  best PIBO config!
seed <- 4

res <- readRDS(paste0("RESULTS_PIBO/results_30_",kern,"_",xi,"_11_500_1.RDS"))
res <- res[res$seed==seed,]
pibo.ix <- which.min(res$y)

X <- matrix(as.numeric(res[pibo.ix, grep("x.",names(res),fixed=T)]),m,2,byrow=T)
P <- P_ + X

if( m==5 ) {
  clrs <- c("deepskyblue", "red", "green3", "purple", "orange2" )
} else {
  clrs <- rainbow(m)
}

par(mar=c(2.1,3.1,2.1,1.1))

plot( P_, pch=21, cex=2, bg=clrs, xlim=c(-1,1), ylim=c(-1,1), xlab="", ylab="", main="PIBO",
      cex.axis=2, cex.main=2 )
rect( 0,0,1,1, border="grey", lwd=2 )
for( i in 1:m )
  lines( c(P_[i,1],P[i,1]), c(P_[i,2],P[i,2]), col=clrs[i], lwd=2 )
legend("topleft",legend=paste("f(P) =",-round(res$y[pibo.ix],2)), cex=2, bty="n" )
points( P, pch=21, cex=2, bg=clrs )


# vanilla BO on flows

res <- readRDS(paste0("RESULTS_vanilla_BO_on_flows/results_30_",kern,"_",xi,"_11_500_1.RDS"))
res <- res[res$seed==seed,]
bofl.ix <- which.min(res$y)

X <- matrix(as.numeric(res[bofl.ix, grep("x.",names(res),fixed=T)]),m,2,byrow=T)
P <- P_ + X

plot( P_, pch=21, cex=2, bg=clrs, xlim=c(-1,1), ylim=c(-1,1), xlab="", ylab="", main="vanilla BO on flows",
      cex.main=2, cex.axis=2 )
rect( 0,0,1,1, border="grey", lwd=2 )
for( i in 1:m )
  lines( c(P_[i,1],P[i,1]), c(P_[i,2],P[i,2]), col=clrs[i], lwd=2 )
legend("topleft",legend=paste("f(P) =",-round(res$y[bofl.ix],2)), cex=2, bty="n" )
points( P, pch=21, cex=2, bg=clrs )



# vanilla BO on point-clouds

res <- readRDS(paste0("RESULTS_vanilla_BO_on_pointclouds/results_30_",kern,"_",xi,"_11_500_1.RDS"))
res <- res[res$seed==seed,]
bopc.ix <- which.min(res$y)

# contrary to PIBO and vanilla-BO-on-flows, vanilla-BO-on-pointclouds works directly in the Physical Space
P <- matrix(as.numeric(res[bopc.ix, grep("x.",names(res),fixed=T)]),m,2,byrow=T)

plot( NA,NA, pch=21, cex=2, bg=clrs, xlim=c(-1,1), ylim=c(-1,1), xlab="", ylab="", main="vanilla BO on point-clouds",
      cex.main=2, cex.axis=2 )
rect( 0,0,1,1, border="grey", lwd=2 )
legend("topleft",legend=paste("f(P) =",-round(res$y[bopc.ix],2)), cex=2, bty="n" )
points( P, pch=21, cex=2, bg=clrs )



#*******************************************************************************************


# PIBO

res <- readRDS(paste0("RESULTS_PIBO/results_30_",kern,"_",xi,"_11_500_1.RDS"))
minixs <- aggregate( res$y, by=list(res$seed), which.min ); names(minixs) <- c("seed","ix")

pibo.ixs <- aggregate(res$y,by=list(res$seed),which.min)
colixs <- grep("x.",names(res),fixed=T)
Xs <- NULL
for( i in 1:nrow(pibo.ixs) )
  Xs <- rbind( Xs, as.numeric(res[pibo.ixs$x[i],colixs]) )

plot( NA, NA, xlim=0:1, ylim=0:1, xlab="", ylab="", main="PIBO", cex.main=2, cex.axis=2 )
for( i in 1:nrow(Xs) )
  points( P_ + matrix(Xs[i,],m,2,byrow=T), pch=21, cex=3, bg=adjustcolor(clrs,alpha.f=0.6) ) 



# vanilla BO on flows

res <- readRDS(paste0("RESULTS_vanilla_BO_on_flows/results_30_",kern,"_",xi,"_11_500_1.RDS"))
minixs <- aggregate( res$y, by=list(res$seed), which.min ); names(minixs) <- c("seed","ix")

bofl.ixs <- aggregate(res$y,by=list(res$seed),which.min)
colixs <- grep("x.",names(res),fixed=T)
Xs <- NULL
for( i in 1:nrow(bofl.ixs) )
  Xs <- rbind( Xs, as.numeric(res[bofl.ixs$x[i],colixs]) )

plot( NA, NA, xlim=0:1, ylim=0:1, main="vanilla BO on flows", cex.main=2, cex.axis=2 )
for( i in 1:nrow(Xs) )
  points( P_ + matrix(Xs[i,],m,2,byrow=T), pch=21, cex=3, bg=adjustcolor(clrs,alpha.f=0.6) ) 



# vanilla BO on point clouds

res <- readRDS(paste0("RESULTS_vanilla_BO_on_pointclouds/results_30_",kern,"_",xi,"_11_500_1.RDS"))
minixs <- aggregate( res$y, by=list(res$seed), which.min ); names(minixs) <- c("seed","ix")

bopc.ixs <- aggregate(res$y,by=list(res$seed),which.min)
colixs <- grep("x.",names(res),fixed=T)
Ps <- NULL
for( i in 1:nrow(bopc.ixs) )
  Ps <- rbind( Ps, as.numeric(res[bopc.ixs$x[i],colixs]) )

plot( NA, NA, xlim=0:1, ylim=0:1, main="vanilla BO on point-clouds", cex.main=2, cex.axis=2 )
for( i in 1:nrow(Ps) )
  points( matrix(Ps[i,],m,2,byrow=T), pch=21, cex=3, bg=adjustcolor(clrs,alpha.f=0.6) ) 
