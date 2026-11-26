rm(list=ls()); graphics.off(); cat("\014")

# 1000 samples for optimizing acquisition
PIBO1000 <- readRDS("additional_RESULTS_PIBO/results_30_exp_6_11_500_1_1000.RDS")
PIBO1000bs <- aggregate(PIBO1000$y,by=list(PIBO1000$seed),min)
PIBO1000tt <- aggregate(PIBO1000$trnTime,by=list(PIBO1000$seed),sum,na.rm=T)
PIBO1000at <- aggregate(PIBO1000$acqTime,by=list(PIBO1000$seed),sum,na.rm=T)
PIBO1000tot <- aggregate(PIBO1000$acqTime+PIBO1000$trnTime,by=list(PIBO1000$seed),sum,na.rm=T)

# 2500 samples for optimizing acquisition
PIBO2500 <- readRDS("additional_RESULTS_PIBO/results_30_exp_6_11_500_1_2500.RDS")
PIBO2500bs <- aggregate(PIBO2500$y,by=list(PIBO2500$seed),min)
PIBO2500tt <- aggregate(PIBO2500$trnTime,by=list(PIBO2500$seed),sum,na.rm=T)
PIBO2500at <- aggregate(PIBO2500$acqTime,by=list(PIBO2500$seed),sum,na.rm=T)
PIBO2500tot <- aggregate(PIBO2500$acqTime+PIBO2500$trnTime,by=list(PIBO2500$seed),sum,na.rm=T)

# 5000 samples for optimizing acquisition
PIBO5000 <- readRDS("additional_RESULTS_PIBO/results_30_exp_6_11_500_1_5000.RDS")
PIBO5000bs <- aggregate(PIBO5000$y,by=list(PIBO5000$seed),min)
PIBO5000tt <- aggregate(PIBO5000$trnTime,by=list(PIBO5000$seed),sum,na.rm=T)
PIBO5000at <- aggregate(PIBO5000$acqTime,by=list(PIBO5000$seed),sum,na.rm=T)
PIBO5000tot <- aggregate(PIBO5000$acqTime+PIBO5000$trnTime,by=list(PIBO5000$seed),sum,na.rm=T)

# (original) 10000 samples for optimizing acquisition
PIBO10k <- readRDS("RESULTS_PIBO/results_30_exp_6_11_500_1.RDS")
PIBO10kbs <- aggregate(PIBO10k$y,by=list(PIBO10k$seed),min)
PIBO10ktt <- aggregate(PIBO10k$trnTime,by=list(PIBO10k$seed),sum,na.rm=T)
PIBO10kat <- aggregate(PIBO10k$acqTime,by=list(PIBO10k$seed),sum,na.rm=T)
PIBO10ktot <- aggregate(PIBO10k$acqTime+PIBO10k$trnTime,by=list(PIBO10k$seed),sum,na.rm=T)


cat("Setting\t\t Obj avg(sd)\t\t trn time avg(sd)\t acq time avg(sd)\t total time avg(sd)\n")
cat("PIBO1000\t",round(mean(PIBO1000bs$x),2),"(",round(sd(PIBO1000bs$x),2),")\t",
    round(mean(PIBO1000tt$x),2),"(",round(sd(PIBO1000tt$x),2),")\t",
    round(mean(PIBO1000at$x),2),"(",round(sd(PIBO1000at$x),2),")\t\t",
    round(mean(PIBO1000tot$x),2),"(",round(sd(PIBO1000tot$x),2),")\n" )
cat("PIBO2500\t",round(mean(PIBO2500bs$x),2),"(",round(sd(PIBO2500bs$x),2),")\t",
    round(mean(PIBO2500tt$x),2),"(",round(sd(PIBO2500tt$x),2),")\t",
    round(mean(PIBO2500at$x),2),"(",round(sd(PIBO2500at$x),2),")\t",
    round(mean(PIBO2500tot$x),2),"(",round(sd(PIBO2500tot$x),2),")\n" )
cat("PIBO5000\t",round(mean(PIBO5000bs$x),2),"(",round(sd(PIBO5000bs$x),2),")\t",
    round(mean(PIBO5000tt$x),2),"(",round(sd(PIBO5000tt$x),2),")\t",
    round(mean(PIBO5000at$x),2),"(",round(sd(PIBO5000at$x),2),")\t",
    round(mean(PIBO5000tot$x),2),"(",round(sd(PIBO5000tot$x),2),")\n" )
cat("PIBO10k\t\t",round(mean(PIBO10kbs$x),2),"(",round(sd(PIBO10kbs$x),2),")\t",
    round(mean(PIBO10ktt$x),2),"(",round(sd(PIBO10ktt$x),2),")\t",
    round(mean(PIBO10kat$x),2),"(",round(sd(PIBO10kat$x),2),")\t",
    round(mean(PIBO10ktot$x),2),"(",round(sd(PIBO10ktot$x),2),")" )


PIBO1000bs <- -PIBO1000bs; PIBO2500bs <- -PIBO2500bs; PIBO5000bs <- -PIBO5000bs; PIBO10kbs <- -PIBO10kbs;

PIBO1000_lo <- mean(PIBO1000bs$x) - sd(PIBO1000bs$x); PIBO1000_up <- mean(PIBO1000bs$x) + sd(PIBO1000bs$x)
PIBO2500_lo <- mean(PIBO2500bs$x) - sd(PIBO2500bs$x); PIBO2500_up <- mean(PIBO2500bs$x) + sd(PIBO2500bs$x)
PIBO5000_lo <- mean(PIBO5000bs$x) - sd(PIBO5000bs$x); PIBO5000_up <- mean(PIBO5000bs$x) + sd(PIBO5000bs$x)
PIBO10k_lo <- mean(PIBO10kbs$x) - sd(PIBO10kbs$x); PIBO10k_up <- mean(PIBO10kbs$x) + sd(PIBO10kbs$x)


par(mar=c(9.1,4.6,1.1,2.1))
plot( c(mean(PIBO1000tot$x),mean(PIBO2500tot$x),mean(PIBO5000tot$x),mean(PIBO10ktot$x)),
      c(mean(PIBO1000bs$x),mean(PIBO2500bs$x),mean(PIBO5000bs$x),mean(PIBO10kbs$x)),
      type="o", pch=19, cex=1.5, lwd=2, col="grey", lty=2, ylab="Generated Power", xlab="",
      cex.lab=2, cex.axis=2, xaxt="n",
      #ylim=range(PIBO1000_lo,PIBO2500_lo,PIBO5000_lo,PIBO10k_lo,PIBO1000_up,PIBO2500_up,PIBO5000_up,PIBO10k_up) )
      ylim=c(70.5,71.5) )
axis(1,at=c(mean(PIBO1000tot$x),mean(PIBO2500tot$x),mean(PIBO5000tot$x),mean(PIBO10ktot$x)),cex.axis=2,las=2)
title(xlab="Time [secs]",line=6,cex.lab=2)

lines( rep(mean(PIBO1000tot$x),2), c(PIBO1000_lo,PIBO1000_up), col="red", lwd=3 )
lines( mean(PIBO1000tot$x)+c(-50,50), rep(PIBO1000_lo,2), col="red", lwd=3 ); lines( mean(PIBO1000tot$x)+c(-50,50), rep(PIBO1000_up,2), col="red", lwd=3 )
lines( rep(mean(PIBO2500tot$x),2), c(PIBO2500_lo,PIBO2500_up), col="green3", lwd=3 )
lines( mean(PIBO2500tot$x)+c(-50,50), rep(PIBO2500_lo,2), col="green3", lwd=3 ); lines( mean(PIBO2500tot$x)+c(-50,50), rep(PIBO2500_up,2), col="green3", lwd=3 )
lines( rep(mean(PIBO5000tot$x),2), c(PIBO5000_lo,PIBO5000_up), col="deepskyblue", lwd=3 )
lines( mean(PIBO5000tot$x)+c(-50,50), rep(PIBO5000_lo,2), col="deepskyblue", lwd=3 ); lines( mean(PIBO5000tot$x)+c(-50,50), rep(PIBO5000_up,2), col="deepskyblue", lwd=3 )
lines( rep(mean(PIBO10ktot$x),2), c(PIBO10k_lo,PIBO10k_up), col="blue", lwd=3 )
lines( mean(PIBO10ktot$x)+c(-50,50), rep(PIBO10k_lo,2), col="blue", lwd=3 ); lines( mean(PIBO10ktot$x)+c(-50,50), rep(PIBO10k_up,2), col="blue", lwd=3 )

points(  c(mean(PIBO1000tot$x),mean(PIBO2500tot$x),mean(PIBO5000tot$x),mean(PIBO10ktot$x)),
         c(mean(PIBO1000bs$x),mean(PIBO2500bs$x),mean(PIBO5000bs$x),mean(PIBO10kbs$x)), pch=21,
         col=c("firebrick","darkgreen","blue","black"), bg=c("red","green3","deepskyblue","blue"), cex=2, lwd=2 )

legend("bottomright",legend=c(1000,2500,5000,10000),col=c("red2","green3","deepskyblue","blue"), pch=19, cex=1.5,
       title="samples for GP-CB optimization" )






par(mar=c(9.1,4.6,1.1,2.1))
plot( round(log10(c(mean(PIBO1000tot$x),mean(PIBO2500tot$x),mean(PIBO5000tot$x),mean(PIBO10ktot$x))),2),
      c(mean(PIBO1000bs$x),mean(PIBO2500bs$x),mean(PIBO5000bs$x),mean(PIBO10kbs$x)),
      type="o", pch=19, cex=1.5, lwd=2, col="grey", lty=2, ylab="Generated Power", xlab="",
      cex.lab=2, cex.axis=2, xaxt="n",
      #ylim=range(PIBO1000_lo,PIBO2500_lo,PIBO5000_lo,PIBO10k_lo,PIBO1000_up,PIBO2500_up,PIBO5000_up,PIBO10k_up) )
      ylim=c(70.5,71.5) )
axis(1,at=round(log10(c(mean(PIBO1000tot$x),mean(PIBO2500tot$x),mean(PIBO5000tot$x),mean(PIBO10ktot$x))),2),cex.axis=2,las=2)
title(xlab="Time [secs] - log10",line=6,cex.lab=2)

lines( rep(round(log10(mean(PIBO1000tot$x)),2),2), c(PIBO1000_lo,PIBO1000_up), col="red", lwd=3 )
lines( round(log10(mean(PIBO1000tot$x)),2)+c(-0.01,0.01), rep(PIBO1000_lo,2), col="red", lwd=3 )
lines( round(log10(mean(PIBO1000tot$x)),2)+c(-0.01,0.01), rep(PIBO1000_up,2), col="red", lwd=3 )

lines( rep(round(log10(mean(PIBO2500tot$x)),2),2), c(PIBO2500_lo,PIBO2500_up), col="green3", lwd=3 )
lines( round(log10(mean(PIBO2500tot$x)),2)+c(-0.01,0.01), rep(PIBO2500_lo,2), col="green3", lwd=3 )
lines( round(log10(mean(PIBO2500tot$x)),2)+c(-0.01,0.01), rep(PIBO2500_up,2), col="green3", lwd=3 )

lines( rep(round(log10(mean(PIBO5000tot$x)),2),2), c(PIBO5000_lo,PIBO5000_up), col="deepskyblue", lwd=3 )
lines( round(log10(mean(PIBO5000tot$x)),2)+c(-0.01,0.01), rep(PIBO5000_lo,2), col="deepskyblue", lwd=3 )
lines( round(log10(mean(PIBO5000tot$x)),2)+c(-0.01,0.01), rep(PIBO5000_up,2), col="deepskyblue", lwd=3 )

lines( rep(round(log10(mean(PIBO10ktot$x)),2),2), c(PIBO10k_lo,PIBO10k_up), col="blue", lwd=3 )
lines( round(log10(mean(PIBO10ktot$x)),2)+c(-0.01,0.01), rep(PIBO10k_lo,2), col="blue", lwd=3 )
lines( round(log10(mean(PIBO10ktot$x)),2)+c(-0.01,0.01), rep(PIBO10k_up,2), col="blue", lwd=3 )

points(  round(log10(c(mean(PIBO1000tot$x),mean(PIBO2500tot$x),mean(PIBO5000tot$x),mean(PIBO10ktot$x))),2),
         c(mean(PIBO1000bs$x),mean(PIBO2500bs$x),mean(PIBO5000bs$x),mean(PIBO10kbs$x)), pch=21,
         col=c("firebrick","darkgreen","blue","black"), bg=c("red","green3","deepskyblue","blue"), cex=2, lwd=2 )

legend("bottomright",legend=c(1000,2500,5000,10000),col=c("red2","green3","deepskyblue","blue"), pch=19, cex=1.5,
       title="samples for GP-CB optimization" )







par(mfrow=c(1,2)); par(mar=c(4.1,4.6,1.1,1.1))
plot( c(1000,2500,5000,10000), c(mean(PIBO1000bs$x),mean(PIBO2500bs$x),mean(PIBO5000bs$x),mean(PIBO10kbs$x)),
      type="l", lwd=2, lty=2, col="grey", ylab="Generated Power", xlab="samples to optimize GP-CB",
      ylim=c(70.5,71.5), cex.axis=1.5, cex.lab=1.5, xaxt="n" )
axis(1,c(1000,2500,5000,10000),labels=c("1k","2.5k","5k","10k"),cex.axis=1.5)

lines( rep(1000,2), c(PIBO1000_lo,PIBO1000_up), col="deepskyblue", lwd=3 )
lines( 1000+c(-100,100), rep(PIBO1000_lo,2), col="deepskyblue", lwd=3 ); lines( 1000+c(-100,100), rep(PIBO1000_up,2), col="deepskyblue", lwd=3 )
lines( rep(2500,2), c(PIBO2500_lo,PIBO2500_up), col="deepskyblue", lwd=3 )
lines( 2500+c(-100,100), rep(PIBO2500_lo,2), col="deepskyblue", lwd=3 ); lines( 2500+c(-100,100), rep(PIBO2500_up,2), col="deepskyblue", lwd=3 )
lines( rep(5000,2), c(PIBO5000_lo,PIBO5000_up), col="deepskyblue", lwd=3 )
lines( 5000+c(-100,100), rep(PIBO5000_lo,2), col="deepskyblue", lwd=3 ); lines( 5000+c(-100,100), rep(PIBO5000_up,2), col="deepskyblue", lwd=3 )
lines( rep(10000,2), c(PIBO10k_lo,PIBO10k_up), col="deepskyblue", lwd=3 )
lines( 10000+c(-100,100), rep(PIBO10k_lo,2), col="deepskyblue", lwd=3 ); lines( 10000+c(-100,100), rep(PIBO10k_up,2), col="deepskyblue", lwd=3 )


points(  c(1000,2500,5000,10000), c(mean(PIBO1000bs$x),mean(PIBO2500bs$x),mean(PIBO5000bs$x),mean(PIBO10kbs$x)), pch=21,
         col="blue", bg="deepskyblue", cex=1.5, lwd=2 )




plot( c(1000,2500,5000,10000), c(mean(PIBO1000tot$x),mean(PIBO2500tot$x),mean(PIBO5000tot$x),mean(PIBO10ktot$x)),
      type="l", lty=2, lwd=2, col="lightgreen", ylab="Time [secs]", xlab="samples to optimize GP-CB",
      ylim=c(100,3500), cex.axis=1.5, cex.lab=1.5, xaxt="n" )
axis(1,c(1000,2500,5000,10000),labels=c("1k","2.5k","5k","10k"),cex.axis=1.5)

lines( rep(1000,2), mean(PIBO1000tot$x)+c(-sd(PIBO1000tot$x),sd(PIBO1000tot$x)), col="green3", lwd=3 )
lines( 1000+c(-100,100), rep(mean(PIBO1000tot$x)-sd(PIBO1000tot$x),2), col="green3", lwd=3 )
lines( 1000+c(-100,100), rep(mean(PIBO1000tot$x)+sd(PIBO1000tot$x),2), col="green3", lwd=3 )

lines( rep(2500,2), mean(PIBO2500tot$x)+c(-sd(PIBO2500tot$x),sd(PIBO2500tot$x)), col="green3", lwd=3 )
lines( 2500+c(-100,100), rep(mean(PIBO2500tot$x)-sd(PIBO2500tot$x),2), col="green3", lwd=3 )
lines( 2500+c(-100,100), rep(mean(PIBO2500tot$x)+sd(PIBO2500tot$x),2), col="green3", lwd=3 )

lines( rep(5000,2), mean(PIBO5000tot$x)+c(-sd(PIBO5000tot$x),sd(PIBO5000tot$x)), col="green3", lwd=3 )
lines( 5000+c(-100,100), rep(mean(PIBO5000tot$x)-sd(PIBO5000tot$x),2), col="green3", lwd=3 )
lines( 5000+c(-100,100), rep(mean(PIBO5000tot$x)+sd(PIBO5000tot$x),2), col="green3", lwd=3 )

lines( rep(10000,2), mean(PIBO10ktot$x)+c(-sd(PIBO10ktot$x),sd(PIBO10ktot$x)), col="green3", lwd=3 )
lines( 10000+c(-100,100), rep(mean(PIBO10ktot$x)-sd(PIBO10ktot$x),2), col="green3", lwd=3 )
lines( 10000+c(-100,100), rep(mean(PIBO10ktot$x)+sd(PIBO10ktot$x),2), col="green3", lwd=3 )



# lines( c(1000,2500,5000,10000), c(mean(PIBO1000tt$x),mean(PIBO2500tt$x),mean(PIBO5000tt$x),mean(PIBO10ktt$x)),
#       type="l", lty=2, lwd=2, col="orange" )
# lines( rep(1000,2), mean(PIBO1000tt$x)+c(-sd(PIBO1000tt$x),sd(PIBO1000tt$x)), col="orange", lwd=3 )
# lines( 1000+c(-100,100), rep(mean(PIBO1000tt$x)-sd(PIBO1000tt$x),2), col="orange", lwd=3 )
# lines( 1000+c(-100,100), rep(mean(PIBO1000tt$x)+sd(PIBO1000tt$x),2), col="orange", lwd=3 )
# 
# lines( rep(2500,2), mean(PIBO2500tt$x)+c(-sd(PIBO2500tt$x),sd(PIBO2500tt$x)), col="orange", lwd=3 )
# lines( 2500+c(-100,100), rep(mean(PIBO2500tt$x)-sd(PIBO2500tt$x),2), col="orange", lwd=3 )
# lines( 2500+c(-100,100), rep(mean(PIBO2500tt$x)+sd(PIBO2500tt$x),2), col="orange", lwd=3 )
# 
# lines( rep(5000,2), mean(PIBO5000tt$x)+c(-sd(PIBO5000tt$x),sd(PIBO5000tt$x)), col="orange", lwd=3 )
# lines( 5000+c(-100,100), rep(mean(PIBO5000tt$x)-sd(PIBO5000tt$x),2), col="orange", lwd=3 )
# lines( 5000+c(-100,100), rep(mean(PIBO5000tt$x)+sd(PIBO5000tt$x),2), col="orange", lwd=3 )
# 
# lines( rep(10000,2), mean(PIBO10ktt$x)+c(-sd(PIBO10ktt$x),sd(PIBO10ktt$x)), col="orange", lwd=3 )
# lines( 10000+c(-100,100), rep(mean(PIBO10ktt$x)-sd(PIBO10ktt$x),2), col="orange", lwd=3 )
# lines( 10000+c(-100,100), rep(mean(PIBO10ktt$x)+sd(PIBO10ktt$x),2), col="orange", lwd=3 )


lines( c(1000,2500,5000,10000), c(mean(PIBO1000at$x),mean(PIBO2500at$x),mean(PIBO5000at$x),mean(PIBO10kat$x)),
       type="l", lty=2, lwd=2, col="red2" )
lines( rep(1000,2), mean(PIBO1000at$x)+c(-sd(PIBO1000at$x),sd(PIBO1000at$x)), col="red2", lwd=3 )
lines( 1000+c(-100,100), rep(mean(PIBO1000at$x)-sd(PIBO1000at$x),2), col="red2", lwd=3 )
lines( 1000+c(-100,100), rep(mean(PIBO1000at$x)+sd(PIBO1000at$x),2), col="red2", lwd=3 )

lines( rep(2500,2), mean(PIBO2500at$x)+c(-sd(PIBO2500at$x),sd(PIBO2500at$x)), col="red2", lwd=3 )
lines( 2500+c(-100,100), rep(mean(PIBO2500at$x)-sd(PIBO2500at$x),2), col="red2", lwd=3 )
lines( 2500+c(-100,100), rep(mean(PIBO2500at$x)+sd(PIBO2500at$x),2), col="red2", lwd=3 )

lines( rep(5000,2), mean(PIBO5000at$x)+c(-sd(PIBO5000at$x),sd(PIBO5000at$x)), col="red2", lwd=3 )
lines( 5000+c(-100,100), rep(mean(PIBO5000at$x)-sd(PIBO5000at$x),2), col="red2", lwd=3 )
lines( 5000+c(-100,100), rep(mean(PIBO5000at$x)+sd(PIBO5000at$x),2), col="red2", lwd=3 )

lines( rep(10000,2), mean(PIBO10kat$x)+c(-sd(PIBO10kat$x),sd(PIBO10kat$x)), col="red2", lwd=3 )
lines( 10000+c(-100,100), rep(mean(PIBO10kat$x)-sd(PIBO10kat$x),2), col="red2", lwd=3 )
lines( 10000+c(-100,100), rep(mean(PIBO10kat$x)+sd(PIBO10kat$x),2), col="red2", lwd=3 )


points(  c(1000,2500,5000,10000), c(mean(PIBO1000tot$x),mean(PIBO2500tot$x),mean(PIBO5000tot$x),mean(PIBO10ktot$x)), pch=21,
         col="darkgreen", bg="green3", cex=1.5, lwd=2 )
# points(  c(1000,2500,5000,10000), c(mean(PIBO1000tt$x),mean(PIBO2500tt$x),mean(PIBO5000tt$x),mean(PIBO10ktt$x)), pch=21,
#          col="orange4", bg="orange", cex=1.5, lwd=2 )
points(  c(1000,2500,5000,10000), c(mean(PIBO1000at$x),mean(PIBO2500at$x),mean(PIBO5000at$x),mean(PIBO10kat$x)), pch=21,
         col="firebrick", bg="red2", cex=1.5, lwd=2 )

# legend("topleft",legend=c("GP fitting","GP-CB optimization","Total"), col=c("orange","red2","green3"), lwd=3, cex=1.5 )
legend("topleft",legend=c("Overall","GP-CB optimization"), col=c("green3","red2"), lwd=3, cex=1.3 )


      