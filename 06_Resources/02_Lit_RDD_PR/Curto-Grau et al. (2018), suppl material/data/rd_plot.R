setwd( "c:\\Users\\mcurto\\Dropbox\\Paper electoral competition, March 2015\\rev") 


getwd() 


## library(extrafont) 
## font_import()
## loadfonts(device="win") 
## fonts()

library(foreign)


rm(list=ls())
data <- read.dta("graphs\\for_LATEplot_CI_epa.dta")
attach(data)
lw<-0.8


pdf("LATEplot_CI_epa.pdf", family="sans", height=5, width=7)
par(mar=c(5.5,5.5, 1.5, 1.5) )
plot(x=B4, y=B1,  type="l",  xlab="Bandwidth", ylab="Regional transfers",lty=1, col.lab='black', cex.lab =0.8, xlim=c(4,30), yaxt="n", xaxt="n", ylim=c(-50,200)) ##, main="Local linear regression with triangular kernel"
abline(h=0, lty=3)
box(col = 'gray25')
axis(1, at=seq(4,30,2), lwd=0.2,  labels=seq(4,30,2), cex.axis =0.8, col.lab='black', cex.lab = 0.8, tck=-0.02 )
axis(2, at=seq(-50,200,50), lwd=0.2,  labels=seq(-50,200,50), cex.axis = 0.8, col.lab='black', cex.lab = 0.8 , tck=-0.02)
lines(x=B4, y=B2,  type="l", lwd=lw, col="black", lty=2)
lines(x=B4, y=B3,  type="l", lwd=lw, col="black", lty=2)
dev.off()



rm(list=ls())
data <- read.dta("graphs\\for_LATEplot_CI_tri.dta")
attach(data)
lw<-0.8

pdf("LATEplot_CI_tri.pdf", family="sans", height=5, width=7)
par(mar=c(5.5,5.5, 1.5, 1.5) )
plot(x=B4, y=B1,  type="l",  xlab="Bandwidth", ylab="Regional transfers",lty=1, col.lab='black', cex.lab =0.8, xlim=c(4,30), yaxt="n", xaxt="n", ylim=c(-50,200)) ##, main="Local linear regression with triangular kernel"
abline(h=0, lty=3)
box(col = 'gray25')
axis(1, at=seq(4,30,2), lwd=0.2,  labels=seq(4,30,2), cex.axis =0.8, col.lab='black', cex.lab = 0.8, tck=-0.02 )
axis(2, at=seq(-50,200,50), lwd=0.2,  labels=seq(-50,200,50), cex.axis = 0.8 , col.lab='black', cex.lab = 0.8, tck=-0.02)
lines(x=B4, y=B2,  type="l", lwd=lw, col="black", lty=2)
lines(x=B4, y=B3,  type="l", lwd=lw, col="black", lty=2)
dev.off()

