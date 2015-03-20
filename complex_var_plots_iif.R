#Project: Moran modelling RT2
#Creator: Teri Davies
#Date created: 04/08/14
#Date last updated: 05/12/14
#Date last run: 05/12/14
#Script purpose: pooling variance results across replicate runs and creation of informative plots [run locally]

#note: if the population changes from 10,000 will need to alter this code 

#manually change directory to the appropriate folder and copy and paste remaining code

library(lokern)

#define the myresamp function for bootstrap approach later
#McMurry T. L., Politis D. N. Bootstrap confidence intervals in nonparametric regression with built-in bias 
#correction. Statistics & Probability Letters 78, 2463-2469 (2008)

myresamp = function(x,residx,est,neighb)
#This function assumes that x is ordered. It takes the input y vector (residx) and finds the nearest neighb neighbours + itself
#It ensures that the window of neighbours is always within 1 to length(x), so at the edges we just extend the window inwards instead
#of having it centered around the target point.
{
	out = numeric(length(x))
	for(j in 1:length(x)){
		iflag = F
		wc = j
		wmin = wc - neighb/2
		wmax = wc + neighb/2
		if(wmin <=0){
			iflag=T
			incr = 0 - wmin+1
			wmin = wmin + incr
			wc = wc + incr
			wmax = wmax + incr
		}
		if(wmax > length(x)){
			if(iflag)stop("problem: input x too short?")
			decr = wmax-length(x)
			wmax = wmax - decr
			wc = wc - decr
			wmin = wmin - decr
			if(wmin <= 0)stop("problem (2): input x too short?")
		}
		out[j] = sample(residx[wmin:wmax],1) + est[j]
		
	}
	out
}

#define the number of individuals in the population
pop<-10000
tot_runs<-100

############################################which model is this?
outfile1="model?_normalised_empirical.txt"
outfile2="model?_analytical.txt"
outfile3="model?_xvalues.txt"
outfile4="model?_normalised_empirical_lowerlim.txt"
outfile5="model?_normalised_empirical_upperlim.txt"


#define the filenames
for (i in 1:tot_runs){
filename=paste("var_plot_",i,".txt",sep="")
data=matrix(scan(filename,what=numeric(), n=(pop-1)*2, sep="\t",na.strings="NA"),(pop-1),2,byrow=TRUE)
if (i==1){
master=data
}
if (i>1){
master=cbind(master,data[,2])
}
}

#the files from the individual runs are (pop-1) by 2 matrices. Each row represents the frequency of 1 to (pop-1) out of pop, in ascending order,
#as you descend the columns. The first column is the frequency and the second is the variance sum for that mutation frequency.
#Variance sum=NA signify missing data (there was no mutation with this frequency for this run)
#To pool together the different runs, average across the variances for each frequency

sum_master=matrix(NA,(pop-1),2)
#insert the frequencies (whole numbers 1-(pop-1))
sum_master[,1]=master[,1]

master_var_only=master[,-1]

for (i in 1:(pop-1)){
if (sum(is.na(master_var_only[i,]))<tot_runs){
sum_master[i,2]=sum(master_var_only[i,],na.rm=TRUE)/(tot_runs)} #average across when there is at least one non missing value
else {
sum_master[i,2]=0 
}
}

#convert the frequencies to decimal values
sum_master[,1]=sum_master[,1]/pop
message("allele frequencies expressed as decimals")

#calculate the variance in the mean estimates
master.zero=ifelse(is.na(master_var_only),0,master_var_only)
v.est=apply(master.zero,1,var,na.rm=T)
#v.est2=ifelse(is.na(v.est),0,v.est)
x1 = log10(c(1:9999)/10000)
y1 = sqrt(v.est/100) #Computes the 'standard error' - the standard deviation in our estimate of the mean

#because the standard error is rather noisy, use lokern to get a smoother estimate 

fit1 = glkerns(x1,y1,hetero=T,is.rand=F)
var.best = fitted(fit1)
var.best = var.best^2 # lokern needs the variance not the standard deviation

############################################change directory to where you want the plot saved

############################################change plot name
a<-"model11_fig_051214"
b<-".tiff"
plotname<-paste(a, b, sep="")
tiff(file=plotname, units="mm", width=xxx, height=xxx, res=300)

#set up for multiple plots in one window
par(mfrow=c(2,2))

#plot the log_10 transformed frequencies against variance contribution
plot(log10(sum_master[,1]),sum_master[,2],xlab="log10(allele frequency)",ylab="mean variance contribution")
title(outer=FALSE,adj=0,main="a",line=1)

#glkerns regression function 
fit1a=glkerns(log10(sum_master[,1]),sum_master[,2],hetero=F,is.rand=F,sig=var.best)


#use bootstrap to obtain upper and lower bounds of variance contribution
#Here we do bootstrap resampling following the outline algorithm in e.g. 
#McMurray and Politis (2008) and papers cited therein. In our case we 
#modify it to sample within a window around each x because of strong 
#heteroscedasticity.

residx = fit1a$y - fit1a$est
bootnum=2000
bootmat = matrix(ncol=length(fit1a$x),nrow=bootnum)
neighb=50
#So here we do bootnum simulations with resampled residuals and 
#put them in bootmat
 

for(j in 1:bootnum){
	pseudo = myresamp(fit1a$x,residx,fit1a$est,neighb)
	fitx = glkerns(fit1a$x,pseudo,hetero=F,is.rand=F,sig=var.best)
	bootmat[j,] = fitx$est
}
q.975 = apply(bootmat,2,quantile,0.975)
q.025 = apply(bootmat,2,quantile,0.025)

#add lines to the previous plot demonstrating the glkerns smoothing with error bars
lines(fit1a$x.out,fit1a$est,col=2,lwd=2)
lines(fit1a$x.out,q.025,col=2,lwd=2,lty=2)
lines(fit1a$x.out,q.975,col=2,lwd=2,lty=2)

#transform the predicted values using the EW (2010) change of variables y=log_10(x)
tf=fit1a$est*log(10)*10^fit1a$x.out

#plot the variance against frequency following the change of variables on the variance axis
plot(log10(sum_master[,1]),tf,type="l",xlab="log10(allele frequency)",ylab="mean variance contribution - change of variables",col=2,lwd=2)
title(outer=FALSE,adj=0,main="b",line=1)
#also plot the transformed error bars
lines(fit1a$x.out,q.025*log(10)*10^fit1a$x.out,col=2,lwd=2,lty=2)
lines(fit1a$x.out,q.975*log(10)*10^fit1a$x.out,col=2,lwd=2,lty=2)

#normalise the plot by estimating the area under the curve and dividing 
#
############################################REMEMBER TO UPDATE YLIM
#
ylim=c(0,0.8)
lfreq=log10(sum_master[,1])
dens=ifelse(tf>=0,tf,0) # restrict variances less than 0 to be 0
tf_l=q.025*log(10)*10^fit1a$x.out
tf_u=q.975*log(10)*10^fit1a$x.out
dens_l=ifelse(tf_l>=0,tf_l,0)
dens_u=ifelse(tf_u>=0,tf_u,0)
d=diff(lfreq)
e=diff(dens)
area=sum(d*e/2)+sum(d*dens[1:length(d)])
plot(lfreq,dens/area,xlab="log10(allele frequency)",ylab="normalised variance contribution - change of variables",ylim=ylim,col=2, cex=0.4)
#plot(lfreq,dens/area,xlab="log10(allele frequency)",ylab="normalised variance contribution - change of variables")
title(outer=FALSE,adj=0,main="c",line=1)
lines(lfreq,dens_l/area,col=2,lwd=2,lty=2)
lines(lfreq,dens_u/area,col=2,lwd=2,lty=2)


#superimpose the analytical result using re-calculated equation [7] from Eyre-Walker (2010) using integration limits 0.0001 and 0.9999

############################################for the original model:

############################################define the scalars
tau=1 # model dependent
sbar=3 # model dependent

beta=0.2
size=10000
vxvt=c(1:(size-1))
x=0
xaxis=c(1:(size-1))
a=0.0001
b=0.9999

#define zeta function
hurwitz=function(s,a,n){
nvec = c(0:n)
res = 1/(a+nvec)^s
sum(res)
}

#test the function - expect 159.524
hurwitz(2.2,0.1 + 0.2/3000,10000)

for (i in 1:(size-1)){
 
x=1/size+x
xaxis[i]=x

#y=log(x, base=10)
 
vxvt_num=hurwitz(2*tau+beta, x+beta/sbar,100000)-hurwitz(2*tau+beta,(sbar+beta)/sbar,100000)
 
vxvt_denom=(1/(2*tau+beta-1))*(hurwitz(2*tau+beta-1,a+beta/sbar,100000)-hurwitz(2*tau+beta-1,b+beta/sbar,100000))+(a-b)*hurwitz(2*tau+beta,(sbar+beta)/sbar,100000)
 
vxvt[i]=2.303*x*vxvt_num/vxvt_denom
}

new_xaxis=log(xaxis,base=10)

#use trapezoidal integration to re-normalise the analytical solution
mdiff=diff(new_xaxis)
adiff=diff(vxvt)
area_analyt=sum(mdiff*adiff/2)+sum(mdiff*vxvt[1:length(mdiff)])
vxvt=vxvt/area_analyt

#add lines to the previous plot showing the analytical result with alteration of integration limits
par(new=T)
plot(new_xaxis,vxvt,type="l",col="green",ylab="",xlab="",ylim=ylim,lwd=2)



############################################for the new model use MB's numerical results:
############################################change model letter of data read in f=mod7, g=mod8, h=mod 9
NUM8=matrix(scan("xxx",what=numeric(), n=9999*3, sep=" "),9999,3,byrow=TRUE)

############################################CHECK NUM8[,1](the log10(allele freqs) equals fit1a$x.out)
#NUM8[,1]==fit1a$x.out - they aren't but random checks show that this is a rounding issue

TFNUM8=NUM8[,2]*log(10)*10^fit1a$x.out
d8=diff(fit1a$x.out)
e8=diff(TFNUM8)
area8=sum(d8*e8/2)+sum(d8*TFNUM8[1:length(d8)])
TFDNUM8=TFNUM8/area8
par(new=T)
plot(fit1a$x.out,TFDNUM8,type="l",col="blue",ylab="",xlab="",ylim=ylim,lwd=2)


#QC
plot(fit1a$x.out,residx,ylab="residuals",xlab="log10(allele frequency)")
abline(h=0)
#fitx = glkerns(log10(sum_master[,1]),residx,hetero=T,is.rand=F)
#lines(fitx$x.out,fitx$est,col="grey",lwd=4)
title(outer=FALSE,adj=0,main="d",line=1)
dev.off()

#plot the analytical result against the raw empirical points - need to transform the analytical result back onto the original scale without change of variables
############################################change plot name
a<-"model11_fig_raw_points_051214"
b<-".tiff"
plotname<-paste(a, b, sep="")
tiff(file=plotname, units="mm", width=xxx, height=xxx, res=300)
par(mfrow=c(2,1))
############################################for models a-e
vxvt_rawscale=vxvt/(xaxis)
plot(vxvt_rawscale,sum_master[,2],xlab="theoretical result divided by allele frequency",ylab="mean raw empirical result")
title(outer=FALSE,adj=0,main="a",line=1)
plot(vxvt,xaxis*sum_master[,2],xlab="theoretical result",ylab="mean raw empirical result multiplied by allele frequency")
title(outer=FALSE,adj=0,main="b",line=1)
############################################for models f-h
plot(TFDNUM8/(10^fit1a$x.out),sum_master[,2],xlab="theoretical result divided by allele frequency",ylab="mean raw empirical result")
title(outer=FALSE,adj=0,main="a",line=1)
plot(TFDNUM8,(10^fit1a$x.out)*sum_master[,2],xlab="theoretical result",ylab="mean raw empirical result multiplied by allele frequency")
title(outer=FALSE,adj=0,main="b",line=1)
dev.off()



############################################write final plot data to file for creating a normalised genetic variance figure across models
############################################need to remove the vxvt row for new model
############################################change working directory
setwd("xxx")
writeit=dens/area
write.table(writeit,outfile1,sep="\t",eol="\n",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(vxvt,outfile2,sep="\t",eol="\n",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(lfreq,outfile3,sep="\t",eol="\n",quote=FALSE,row.names=FALSE, col.names=FALSE)

#write lower limit to file
writeit_l=dens_l/area
write.table(writeit_l,outfile4,sep="\t",eol="\n",quote=FALSE,row.names=FALSE, col.names=FALSE)
#write upper limit to file
writeit_u=dens_u/area
write.table(writeit_u,outfile5,sep="\t",eol="\n",quote=FALSE,row.names=FALSE, col.names=FALSE)
