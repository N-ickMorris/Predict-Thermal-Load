#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# MY FUNCTIONS

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# A function to calculate various transformations on a given 'vector' and ranks transformations by highest magnitude of correlation with the given 'Response' vector

Transform=function(Response,vector)
{
  Results=matrix(ncol=1)
  name=deparse(substitute(vector))

  Results=rbind(Results,abs(cor(logb(vector,exp(1)),Response)),abs(cor(logb(vector,2),Response)),abs(cor(Response,logb(vector,3))),abs(cor(Response,logb(vector,4))),abs(cor(Response,logb(vector,6))),abs(cor(Response,logb(vector,8))),abs(cor(Response,logb(vector,10))),abs(cor(Response,cos(vector))),abs(cor(Response,sin(vector))),abs(cor(Response,tan(vector))),abs(cor(Response,acos(vector))),abs(cor(Response,asin(vector))),abs(cor(Response,atan(vector))),abs(cor(Response,exp(vector))),abs(cor(Response,10^(vector))),abs(cor(Response,vector^-4)),abs(cor(Response,vector^-3)),abs(cor(Response,vector^-2.5)),abs(cor(Response,vector^-2.4)),abs(cor(Response,vector^-2.3)),abs(cor(Response,vector^-2.2)),abs(cor(Response,vector^-2.1)),abs(cor(Response,vector^-2)),abs(cor(Response,vector^-1.9)),abs(cor(Response,vector^-1.8)),abs(cor(Response,vector^-1.7)),abs(cor(Response,vector^-1.6)),abs(cor(Response,vector^-1.5)),abs(cor(Response,vector^-1.4)),abs(cor(Response,vector^-1.3)),abs(cor(Response,vector^-1.2)),abs(cor(Response,vector^-1.1)),abs(cor(Response,vector^-1)),abs(cor(Response,vector^-0.9)),abs(cor(Response,vector^-0.8)),abs(cor(Response,vector^-0.7)),abs(cor(Response,vector^-0.6)),abs(cor(Response,vector^-0.5)),abs(cor(Response,vector^-0.4)),abs(cor(Response,vector^-0.3)),abs(cor(Response,vector^-0.2)),abs(cor(Response,vector^-0.1)),abs(cor(Response,vector^0.1)),abs(cor(Response,vector^0.2)),abs(cor(Response,vector^0.3)),abs(cor(Response,vector^0.4)),abs(cor(Response,vector^0.5)),abs(cor(Response,vector^0.6)),abs(cor(Response,vector^0.7)),abs(cor(Response,vector^0.8)),abs(cor(Response,vector^0.9)),abs(cor(Response,vector^1)),abs(cor(Response,vector^1.1)),abs(cor(Response,vector^1.2)),abs(cor(Response,vector^1.3)),abs(cor(Response,vector^1.4)),abs(cor(Response,vector^1.5)),abs(cor(Response,vector^1.6)),abs(cor(Response,vector^1.7)),abs(cor(Response,vector^1.8)),abs(cor(Response,vector^1.9)),abs(cor(Response,vector^2)),abs(cor(Response,vector^2.1)),abs(cor(Response,vector^2.2)),abs(cor(Response,vector^2.3)),abs(cor(Response,vector^2.4)),abs(cor(Response,vector^2.5)),abs(cor(Response,vector^3)),abs(cor(Response,vector^4)))

  Results=as.matrix(Results[-1,])

  rownames(Results)=c(paste('logb(',name,')',sep=''),paste('logb(',name,',2)',sep=''),paste('logb(',name,',3)',sep=''),paste('logb(',name,',4)',sep=''),paste('logb(',name,',6)',sep=''),paste('logb(',name,',8)',sep=''),paste('logb(',name,',10)',sep=''),paste('cos(',name,')',sep=''),paste('sin(',name,')',sep=''),paste('tan(',name,')',sep=''),paste('acos(',name,')',sep=''),paste('asin(',name,')',sep=''),paste('atan(',name,')',sep=''),paste('exp(',name,')',sep=''),paste('10^',name,sep=''),paste(name,'^',-4,sep=''),paste(name,'^',-3,sep=''),paste(name,'^',-2.5,sep=''),paste(name,'^',-2.4,sep=''),paste(name,'^',-2.3,sep=''),paste(name,'^',-2.2,sep=''),paste(name,'^',-2.1,sep=''),paste(name,'^',-2,sep=''),paste(name,'^',-1.9,sep=''),paste(name,'^',-1.8,sep=''),paste(name,'^',-1.7,sep=''),paste(name,'^',-1.6,sep=''),paste(name,'^',-1.5,sep=''),paste(name,'^',-1.4,sep=''),paste(name,'^',-1.3,sep=''),paste(name,'^',-1.2,sep=''),paste(name,'^',-1.1,sep=''),paste(name,'^',-1,sep=''),paste(name,'^',-0.9,sep=''),paste(name,'^',-0.8,sep=''),paste(name,'^',-0.7,sep=''),paste(name,'^',-0.6,sep=''),paste(name,'^',-0.5,sep=''),paste(name,'^',-0.4,sep=''),paste(name,'^',-0.3,sep=''),paste(name,'^',-0.2,sep=''),paste(name,'^',-0.1,sep=''),paste(name,'^',0.1,sep=''),paste(name,'^',0.2,sep=''),paste(name,'^',0.3,sep=''),paste(name,'^',0.4,sep=''),paste(name,'^',0.5,sep=''),paste(name,'^',0.6,sep=''),paste(name,'^',0.7,sep=''),paste(name,'^',0.8,sep=''),paste(name,'^',0.9,sep=''),paste(name,'^',1,sep=''),paste(name,'^',1.1,sep=''),paste(name,'^',1.2,sep=''),paste(name,'^',1.3,sep=''),paste(name,'^',1.4,sep=''),paste(name,'^',1.5,sep=''),paste(name,'^',1.6,sep=''),paste(name,'^',1.7,sep=''),paste(name,'^',1.8,sep=''),paste(name,'^',1.9,sep=''),paste(name,'^',2,sep=''),paste(name,'^',2.1,sep=''),paste(name,'^',2.2,sep=''),paste(name,'^',2.3,sep=''),paste(name,'^',2.4,sep=''),paste(name,'^',2.5,sep=''),paste(name,'^',3,sep=''),paste(name,'^',4,sep=''))

  Results=Results[order(Results[,1],decreasing = TRUE),]
  Results=t(t(Results))
  colnames(Results)=c("Correlation.Magnitude")
  return(Results)
}

#------------------

# A function to calculate R^2, R^2-Adj, R^2-Pred, and the Absolute Difference between R^2 and R^2-Adj

Rsq.Sum=function(Model)
{
  #Calculate the Predictive Residuals of 'model'
  PR=residuals(Model)/(1-lm.influence(Model)$hat)

  #Calculate the Predicted Residual Sum of Squares
  PRESS=sum(PR^2)

  #Calculate the Total Sum of Squares
  TSS=sum(anova(Model)$"Sum Sq")

  #R^2 Summary
  return(c("R^2"=round(summary(Model)$r.squared,4),"R^2-Adj"=round(summary(Model)$adj.r.squared,4),"R^2-Pred"=round(1-(PRESS/TSS),4),"Diff"=abs(round(summary(Model)$r.squared,4)-round(summary(Model)$adj.r.squared,4))))
}

#------------------

# A function to create 6 residual plots for evaluating the 6 major model assumptions

ResPlots=function(LM)
{
  layout(mat = matrix(1:6,2,3, byrow=FALSE),height = c(3,3))
  par(mar=c(5.1, 5.1, 5.1, 2.1))
  plot(LM,which=ifelse(length(LM$residuals)>30,1,3),pch=16,cex=.75,col="gray27",lwd=2,panel.first=c(rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col="lightgrey"),grid(col="white",lty=1)))
  abline(h=0,col="black",lty=2)
  plot(LM,which=2,pch=16,cex=.75,col="gray27",lwd=2,panel.first=c(rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col="lightgrey"),grid(col="white",lty=1)))
  plot(LM$residuals[-1],LM$residuals[-length(LM$residuals)],font.main=1,cex.main=1.5,main="Residuals vs Residuals - Lag 1",xlab="Residuals without First Observation",ylab="Residuals without Last Observation",pch=16,cex=.75,col="gray27",lwd=2,panel.first=c(rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col="lightgrey"),grid(col="white",lty=1)))
  lines(loess.smooth(LM$residuals[-1],LM$residuals[-length(LM$residuals)]),col="red",lwd=2)
  plot(LM,which=5,pch=16,cex=.75,col="gray27",lwd=2,panel.first=c(rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col="lightgrey"),grid(col="white",lty=1)))
  abline(h=0,col="black",lty=2)
  abline(v=0,col="black",lty=2)
  z=boxplot(LM$residuals,whisklty=1,staplelty=0,notch=TRUE,horizontal=TRUE,outcol="red",outpch=16,ylim=c(min(LM$residuals),max(LM$residuals)))
  z
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col="lightgrey")
  grid(col="white",lty=1)
  boxplot(LM$residuals,add=TRUE,medcol=ifelse(z$conf[1,]<0&&z$conf[2,]>0,"green","red"),whisklty=1,staplelty=0,notch=TRUE,horizontal=TRUE,outcol="red",outpch=16,outcex=.75,ylim=c(min(LM$residuals),max(LM$residuals)))
  title(font.main=1,cex.main=1.5,main=paste("Boxplot of Residuals"))
  points(x=mean(LM$residuals),y=1,col="blue",pch=10,cex=1.5,lwd=1.5)
  plot(density(LM$residuals),font.main=1,cex.main=1.5,main="Kernel Density Plot of Residuals",xlab= " ",col="black",xlim=c(min(LM$residuals),max(LM$residuals)))
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col="lightgrey")
  grid(col="white",lty=1)
  polygon(density(LM$residuals),col="white",xlim=c(min(LM$residuals),max(LM$residuals)))
  lines(density(LM$residuals),lty=1,xlim=c(min(LM$residuals),max(LM$residuals)))
}

#------------------

# A function to plot the residuals v. residuals based on a given 'Lag' to evaluate the major model assumption of 'The residuals are non-correlated'

DW.Plot=function(LM,Lag)
{
  plot(LM$residuals[-c(1:Lag)],LM$residuals[-c((length(LM$residuals)-Lag+1):length(LM$residuals))],font.main=1,cex.main=1.5,main=paste("Residuals vs Residuals - Lag",as.character(Lag)),xlab="Residuals without First 'Lag' Observations",ylab="Residuals without Last 'Lag' Observations",pch=16,cex=.75,col="gray27",lwd=2,panel.first=c(rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col="lightgrey"),grid(col="white",lty=1)))
  lines(loess.smooth(LM$residuals[-c(1:Lag)],LM$residuals[-c((length(LM$residuals)-Lag+1):length(LM$residuals))]),col="red",lwd=2)
  print(noquote("Autocorrelation:"))
  return(cor(LM$residuals[-c(1:Lag)],LM$residuals[-c((length(LM$residuals)-Lag+1):length(LM$residuals))]))
}

#------------------

# A function to analytically evaluate each of the 6 major assumptions

AssumTables=function(LM)
{
  tab1=cbind("N:"=length(LM$residuals),"Min:"=round(min(LM$residuals),6),"1st Qu:"=round(quantile(LM$residuals,.25),6),"Median:"=round(median(LM$residuals),6),"Mean:"=round(mean(LM$residuals),6),"3rd Qu:"=round(quantile(LM$residuals,.75),6),"Max:"=round(max(LM$residuals),6),"Range:"=round(max(LM$residuals)-min(LM$residuals),6),"StdDev:"=round(sd(LM$residuals),6),"IQR:"=round(quantile(LM$residuals,.75)-quantile(LM$residuals,.25),6))
  row.names(tab1)=c("")
  tab2=t.test(LM$residuals)
  require(car)
  tab3=ncvTest(LM)
  tab4=shapiro.test(LM$residuals)
  tab5=durbinWatsonTest(LM)
  tab6=summary(LM)
  require(car)
  tab7=vif(LM)
  Results=list("Basic Summary of Residuals:"=tab1,"Assumption 1 - Residuals have an Expected Value of Zero - Zero Within CI"=tab2,"Assumption 2 - Residuals have Constant Variance - P-Value > 0.05"=tab3,"Assumption 3 - Residuals are Normally Distributed - P-Value > 0.05"=tab4,"Assumption 4 - Residuals are Uncorrelated - P-value > 0.05"=tab5,"Assumption 5 - The Relationship between the Response and Regressors is Correct - All P-Value < 0.1"=tab6,"Assumption 6 - The Regressors are Independent - All vif < 10"=tab7)
  return(Results)
}

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# DATA SET ADJUSTMENTS

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

DS=read.csv(file.choose())

# Reposition response variable

DS=data.frame("TL"=DS[,8],DS[,1:7])

# Fix a Roof value from 1470 to 147

DS[9,5]=147

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# CONTINUOUS VARIABLE CORRELATIONS

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Correlations of Compactness Transformations

head(Transform(DS$TL,DS$Compactness))

# Compactness
# 0.1358848

# I(Compactness^-4)
# -0.2113375

# Piecewise Regression with a break point at .90889

#------------------

# Correlations of Surface Transformations

head(Transform(DS$TL,DS$Surface))

# Surface
# -0.1358848

# tan(Surface)
# 0.5423143

# sin(Surface)
# -0.5300436

# piecewise regresion with a break point at 560

#------------------

# Correlations of Wall Transformations

head(Transform(DS$TL,DS$Wall))

# Wall
# 0.09212918

#------------------

# Correlations of Roof Transformations

head(Transform(DS$TL,DS$Roof))

# Roof
# -0.1881579

# cos(Roof)
# -0.38625

# piecewise regression at 120

#------------------

# Correlations of Glazing Transformations

head(Transform(DS$TL,DS$Glazing))

# Glazing
# 0.785362

# exp(Glazing)
# 0.7868553

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# VARIABLE GRAPHICS

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

require(ggplot2)

# Surface:Compactness:Glazing

ggplot(DS,aes(x=Surface,y=TL))+geom_point(aes(color=Compactness),cex=3)+facet_wrap(~Glazing)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")

#------------------

# Roof:Surface:Glazing

ggplot(DS,aes(x=Roof,y=TL))+geom_point(aes(color=Surface),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Glazing)

#------------------

# Surface:Compactness

ggplot(DS,aes(x=Compactness,y=TL))+geom_point(aes(color=Surface),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Surface)

#------------------

# Wall:Compactness

ggplot(DS,aes(x=Compactness,y=TL))+geom_point(aes(color=Wall),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")

#------------------

# Roof:Compactness

ggplot(DS,aes(x=Compactness,y=TL))+geom_point(aes(color=Roof),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Roof)

#------------------

# Glazing:Compactness

ggplot(DS,aes(x=Compactness,y=TL))+geom_point(aes(color=Glazing),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Glazing)

#------------------

# Surface:Wall

ggplot(DS,aes(x=Surface,y=TL))+geom_point(aes(color=Wall),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")

#------------------

# Surface:Roof

ggplot(DS,aes(x=Surface,y=TL))+geom_point(aes(color=Roof),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Roof)

#------------------

# Surface:Glazing

ggplot(DS,aes(x=Surface,y=TL))+geom_point(aes(color=Glazing),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Glazing)

#------------------

# Wall:Roof

ggplot(DS,aes(x=Wall,y=TL))+geom_point(aes(color=Roof),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Roof)

#------------------

# Wall:Glazing

ggplot(DS,aes(x=Wall,y=TL))+geom_point(aes(color=Glazing),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Glazing)

#------------------

# Roof:Glazing

ggplot(DS,aes(x=Roof,y=TL))+geom_point(aes(color=Glazing),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Glazing)

#------------------

# Orientation:Glazing

ggplot(DS,aes(x=Orientation,y=TL))+geom_point(aes(color=Glazing),cex=3)+scale_color_gradient(low = "#0000FF", high = "#FF0000",na.value = "#00FF00", guide = "legend")+facet_wrap(~Glazing)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# MANUAL MODEL STEPPING

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

lm1=lm(TL~.,data=DS)
lm2=lm(TL~.-Surface,data=DS)
lm3=lm(TL~.-Surface-WindowDist-Orientation,data=DS)
lm4=lm(TL~.-Surface-WindowDist-Orientation-Roof,data=DS)
lm4=lm(TL ~ Compactness + Wall + Glazing, data = DS)
lm5=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface, data = DS)
lm6=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing, data = DS)
lm7=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Wall:Glazing, data = DS)
lm8=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Roof:Glazing, data = DS)
lm9=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Roof:Glazing, data = DS)
lm10=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Orientation:Glazing, data = DS)
lm11=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Orientation:Glazing + Glazing:WindowDist, data = DS)
lm12=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Glazing:WindowDist, data = DS)
lm13=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Orientation:Glazing + Compactness:Surface:Orientation, data = DS)
lm14=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Orientation, data = DS)
lm15=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Orientation + Compactness:Surface:Glazing, data = DS)
lm16=lm(TL ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing, data = DS)
lm17=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Orientation + Compactness:Surface:Glazing, data = DS)
lm18=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Orientation + Compactness:Surface:Glazing + Compactness:Surface:WindowDist, data = DS)
lm19=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Compactness:Surface:WindowDist, data = DS)
lm20=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing, data = DS)
lm21=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Compactness:Wall:Glazing, data = DS)
lm22=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Compactness:Roof:Glazing, data = DS)
lm23=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Roof:Glazing, data = DS)
lm24=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Compactness:Orientation:Glazing, data = DS)
lm25=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Compactness:Orientation:Glazing + Compactness:Glazing:WindowDist, data = DS)
lm26=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Compactness:Glazing:WindowDist, data = DS)
lm27=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Compactness:Orientation:Glazing + Surface:Wall:Glazing, data = DS)
lm28=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Surface:Wall:Glazing, data = DS)
lm29=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Compactness:Orientation:Glazing + Surface:Roof:Glazing, data = DS)
lm30=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Surface:Roof:Glazing, data = DS)
lm31=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Surface:Roof:Glazing, data = DS)
lm32=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Compactness:Orientation:Glazing + Surface:Orientation:Glazing, data = DS)
lm33=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Surface:Orientation:Glazing, data = DS)
lm34=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Surface:Glazing:WindowDist, data = DS)
lm35=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Surface:Glazing:WindowDist, data = DS)
lm36=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing, data = DS)
lm37=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Wall:Roof:Glazing, data = DS)
lm38=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing, data = DS)
lm39=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing, data = DS)
lm40=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing, data = DS)
lm41=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Orientation:Glazing, data = DS)
lm42=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing, data = DS)
lm43=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Wall:Orientation:Glazing, data = DS)
lm44=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Wall:Glazing:WindowDist, data = DS)
lm45=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Glazing:WindowDist, data = DS)
lm46=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Wall:Glazing:WindowDist, data = DS)
lm47=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Wall:Glazing:WindowDist, data = DS)
lm48=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Glazing:WindowDist, data = DS)
lm49=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Wall:Glazing:WindowDist, data = DS)
lm50=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Wall:Glazing:WindowDist, data = DS)
lm51=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Roof:Orientation:Glazing, data = DS)
lm52=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Roof:Orientation:Glazing, data = DS)
lm53=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Roof:Orientation:Glazing, data = DS)
lm54=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Roof:Glazing:WindowDist, data = DS)
lm55=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Roof:Glazing:WindowDist, data = DS)
lm56=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Roof:Glazing:WindowDist, data = DS)
lm57=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Roof:Glazing:WindowDist, data = DS)
lm58=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Roof:Glazing:WindowDist, data = DS)
lm59=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Roof:Glazing:WindowDist, data = DS)
lm60=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Roof:Glazing:WindowDist, data = DS)
lm61=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Orientation:Glazing:WindowDist, data = DS)
lm62=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Orientation:Glazing:WindowDist, data = DS)
lm63=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Compactness:Surface:Wall:Orientation, data = DS)
lm64=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Wall:Orientation, data = DS)
lm65=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Wall:Orientation, data = DS)
lm66=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Surface:Roof:Glazing, data = DS)
lm67=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Compactness:Surface:Roof:Glazing, data = DS)
lm68=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Roof:Glazing, data = DS)
lm69=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Surface:Orientation:Glazing, data = DS)
lm70=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Compactness:Surface:Orientation:Glazing, data = DS)
lm71=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Orientation:Glazing, data = DS)
lm72=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Surface:Glazing:WindowDist, data = DS)
lm73=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Compactness:Surface:Glazing:WindowDist, data = DS)
lm74=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Compactness:Surface:Glazing:WindowDist, data = DS)
lm75=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Glazing:WindowDist, data = DS)
lm76=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing, data = DS)
lm77=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing + Wall:Orientation:Glazing + Compactness:Wall:Roof:Glazing, data = DS)
lm78=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing, data = DS)
lm79=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Orientation:Glazing + Compactness:Wall:Roof:Glazing, data = DS)
lm80=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Orientation:Glazing + Compactness:Wall:Roof:Glazing, data = DS)
lm81=lm(TL ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Wall:Roof:Glazing, data = DS)
lm82=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Compactness:Wall:Orientation:Glazing, data = DS)
lm83=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Compactness:Wall:WindowDist:Glazing, data = DS)
lm84=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Compactness:Roof:Orientation:Glazing, data = DS)
lm85=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Compactness:Roof:WindowDist:Glazing, data = DS)
lm86=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Surface:Roof:Wall:Glazing, data = DS)
lm87=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Surface:Roof:Wall:Glazing, data = DS)
lm88=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Surface:Orientation:Wall:Glazing, data = DS)
lm89=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Surface:WindowDist:Wall:Glazing, data = DS)
lm90=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Surface:Orientation:Roof:Glazing, data = DS)
lm91=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Surface:Orientation:Roof:Glazing, data = DS)
lm92=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Surface:WindowDist:Roof:Glazing, data = DS)
lm93=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Surface:WindowDist:Roof:Glazing, data = DS)
lm94=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Wall:Orientation:Roof:Glazing, data = DS)
lm95=lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing + Wall:WindowDist:Roof:Glazing, data = DS)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# MODEL CANIDATES

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Chosen after the manual model stepping

lmday1 = lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing, data = DS)

# chosen after evaluating 7776 X transformation combinations on lmday1

lm3033=lm(TL~I(Compactness^-1)+I(Compactness^-1):I(Surface^2)+I(Compactness^-1):I(Surface^2):Glazing+I(Compactness^-1):I(Surface^2):I(Wall^-1):Orientation+I(Compactness^-1):I(Wall^-1):I(Roof^-1):Glazing,data=DS)

# chosen after doing a PCA linear combination of Wall and Roof: "Z1"=(-.34*DS$Wall) + (-.94*DS$Roof);"Z2"=(.94*DS$Wall) + (-.34*DS$Roof)

DS2=data.frame(DS[,c(1,2,3)],"Z1"=(-.34*DS$Wall) + (-.94*DS$Roof),"Z2"=(.94*DS$Wall) + (-.34*DS$Roof),DS[,c(6,7,8)])
lmtwopca5=lm(TL ~ (. - Orientation - WindowDist - Surface - Glazing + exp(Glazing))^2 + Orientation - Z1:Z2 - Compactness:Z2 - Compactness:exp(Glazing) - Z2:exp(Glazing), data = DS2)

# chosen after doing a mixed stepwise with a scope of threeway interactions and removing terms with drop1 analysis

lmthree.t=lm(log(TL) ~ Compactness + Orientation + Surface:Glazing + Surface:Roof + Glazing:Compactness + Surface:Compactness + Surface:Compactness:Glazing, data = DS)

# chosen to normalize the residuals of lm3033 by log transforming the response

lm3033.t=lm(log(TL)~I(Compactness^-1)+I(Compactness^-1):I(Surface^2)+I(Compactness^-1):I(Surface^2):Glazing+I(Compactness^-1):I(Surface^2):I(Wall^-1):Orientation+I(Compactness^-1):I(Wall^-1):I(Roof^-1):Glazing,data=DS)

# chosen to normalize the residuals of lm65 from manual model steping by log transforming the response

lm65.t=lm(log(TL) ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Compactness:Surface:Wall:Orientation, data = DS)

# chosen to normalize the residuals of lm15 from manual model steping by log transforming the response

lm15.t=lm(log(TL) ~ Compactness + Wall + Glazing + Compactness:Surface + Compactness:Glazing + Compactness:Surface:Glazing,data=DS)

# chosen to normalize the residuals of lm39 from manual model steping by log transforming the response

lm39.t=lm(log(TL) ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Surface:Orientation:Glazing + Wall:Roof:Glazing, data = DS)

# chosen to normalize the residuals of lm60 from manual model steping by log transforming the response

lm60.t=lm(log(TL) ~ Compactness + Wall + Compactness:Surface + Compactness:Surface:Glazing + Wall:Roof:Glazing + Roof:Glazing:WindowDist, data = DS)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# CHOSEN MODEL

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

lmday1 = lm(TL ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing, data = DS)

# use boxcox to find a power within the 95% CI to improve the normality assumption

boxCox(lmday1,grid=FALSE,panel.first=c(rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col="lightgrey"),grid(col="white",lty=1)),las=1)
lmday1.t = lm(sqrt(TL) ~ Compactness + Compactness:Surface + Compactness:Surface:Glazing + Compactness:Surface:Wall:Orientation + Compactness:Wall:Roof:Glazing, data = DS)

# adjust the terms in the lmday1.t model that violate multicollinearity

vif(lmday1.t)
lmday3=lm(sqrt(TL) ~ Compactness + Orientation + Compactness:Surface + Compactness:Surface:Wall + Compactness:Glazing + Roof:Glazing:Wall, data = DS)
vif(lmday3)

# Chosen Model

lmCHOSEN=lm(sqrt(TL) ~ Compactness + Orientation + Compactness:Surface + Compactness:Surface:Wall + Compactness:Glazing + Roof:Glazing:Wall, data = DS)
