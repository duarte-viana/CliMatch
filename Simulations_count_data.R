###########################################
# Supporting Information
# "Ecological traits linking climate with bird abundances and distributions"
# Duarte S. Viana, Jonathan M. Chase
###########################################

# Simulations

library(performance)

N <- 1000
x <- runif(N,0,1)
y0 <- exp((-(0.5-x)^2)/0.1)
plot(x,y0)

r2 <- c()
d2 <- c()
mean_abund <- c()
dec <- seq(2,100,1)
quartz(height=3,width=10)
par(mfrow=c(1,3),mar=c(3,3,1,1),cex.axis=1.3)
newx <- data.frame(x=seq(0,1,0.01))
for(j in dec){
  y <- j*y0
  y <- rpois(N,y)
  fitj <- glm(y~poly(x,2),family=poisson)
  preds <- predict(fitj, newdata=newx, type="response")
  r2 <- c(r2,cor(y,predict(fitj,type="response"),method="spearman")^2)
  d2 <- c(d2,r2_mcfadden(fitj)$R2)
  mean_abund <- c(mean_abund,mean(y))
  if(j %in% c(5,50,100)){
    plot(x,y,pch=16,xlab="Predictor",ylab="Abundance")
    lines(newx$x,preds,lwd=3,col="red")
  }
}

quartz(height=4,width=4)
par(mar=c(4,4,1,1))
plot(log(mean_abund+1),r2,pch=16,ylim=c(0,1),col="blue",ylab="Fitting performance",xlab="Mean abundance (log)")
points(log(mean_abund+1),d2,pch=16,col="darkred")
legend("bottomright",c("Spearman's","McFadden's"),pch=16,col=c("blue","darkred"),bty="n")
