#Project: Moran modelling RT2
#Creator: Teri Davies
#Date created: 04/08/14
#Date last updated: 28/12/14
#Date last run: 28/12/14
#Script purpose: to create a figure detailing avg kurtosis of trait distribution across runs of each model
#				 and variance weighted mean frequency across runs of each model


#Load R moments package
library("moments")

#READ IN THE VARIANCE DENSITY DATA 
setwd("xxx")
#read in the data by model
infile1="model1_normalised_empirical.txt"
infile3="model1_xvalues.txt"
EMP1=scan(infile1,what=numeric(), n=9999)
AF1=scan(infile3,what=numeric(), n=9999)

setwd("xxx")
#read in the data by model
infile1="model2_normalised_empirical.txt"
infile3="model2_xvalues.txt"
EMP2=scan(infile1,what=numeric(), n=9999)
AF2=scan(infile3,what=numeric(), n=9999)

setwd("xxx")
#read in the data by model
infile1="model3_normalised_empirical.txt"
infile3="model3_xvalues.txt"
EMP3=scan(infile1,what=numeric(), n=9999)
AF3=scan(infile3,what=numeric(), n=9999)

setwd("xxx")
#read in the data by model
infile1="model4_normalised_empirical.txt"
infile3="model4_xvalues.txt"
EMP4=scan(infile1,what=numeric(), n=9999)
AF4=scan(infile3,what=numeric(), n=9999)

setwd("xxx")
#read in the data by model
infile1="model5_normalised_empirical.txt"
infile3="model5_xvalues.txt"
EMP5=scan(infile1,what=numeric(), n=9999)
AF5=scan(infile3,what=numeric(), n=9999)

setwd("xxx")
#read in the data by model
infile1="270814_model7_normalised_empirical.txt"
infile3="270814_model7_xvalues.txt"
EMP7=scan(infile1,what=numeric(), n=9999)
AF7=scan(infile3,what=numeric(), n=9999)

setwd("xxx")
#read in the data by model
infile1="270814_model8_normalised_empirical.txt"
infile3="270814_model8_xvalues.txt"
EMP8=scan(infile1,what=numeric(), n=9999)
AF8=scan(infile3,what=numeric(), n=9999)

setwd("xxx")
#read in the data by model
infile1="270814_model9_normalised_empirical.txt"
infile3="270814_model9_xvalues.txt"
EMP9=scan(infile1,what=numeric(), n=9999)
AF9=scan(infile3,what=numeric(), n=9999)



#READ IN TRAIT DISTN DATA

folder="xxx"
setwd(folder)

kurt_by_mod=c(1:8)
for (k in 1:8){
kurt_by_mod[k]=NA
}

for (mod in c(1,2,3,4,5,7,8,9)){
trait_matrix=matrix(NA,nrow=10000,ncol=100)

for (run in c(1:100)){
filename=paste(folder,"model",mod,"/",run,"_run.txt",sep="")
trait_matrix[,run]=scan(file=filename,what="raw",n=10000*1,sep="\t")
}
mode(trait_matrix)="numeric"

kurtosisval=c(1:100)
for (j in 1:100){
kurtosisval[j]=NA
}

for (i in c(1:100)){

kurtosisval[i]=kurtosis(trait_matrix[,i])
}
if (mod<6){
kurt_by_mod[mod]=mean(kurtosisval)
}
else if (mod>6){
kurt_by_mod[mod-1]=mean(kurtosisval)
}



}

#CD TO KURTOSIS FOLDER
setwd(xxx)

#calculate sum(variance*frequency) from the variance density plots of each model
varFreq=c(1:8)
for (i in 1:8){
varFreq[i]=NA
}

VF1=EMP1*AF1
VF2=EMP2*AF2
VF3=EMP3*AF3
VF4=EMP4*AF4
VF5=EMP5*AF5
VF7=EMP7*AF7
VF8=EMP8*AF8
VF9=EMP9*AF9

varFreq[1]=sum(VF1)
varFreq[2]=sum(VF2)
varFreq[3]=sum(VF3)
varFreq[4]=sum(VF4)
varFreq[5]=sum(VF5)
varFreq[6]=sum(VF7)
varFreq[7]=sum(VF8)
varFreq[8]=sum(VF9)

#calculate sum(variance*frequency)/sum(variance) 
varsum=c(1:8)
for (i in 1:8){
varsum[i]=NA
}

varsum[1]=sum(EMP1)
varsum[2]=sum(EMP2)
varsum[3]=sum(EMP3)
varsum[4]=sum(EMP4)
varsum[5]=sum(EMP5)
varsum[6]=sum(EMP7)
varsum[7]=sum(EMP8)
varsum[8]=sum(EMP9)

varFreq_dVar=varFreq/varsum

#create scatter plot with kurtosis on the y-axis, filled points for EW model and hollow points for new model
tiff(file="mean_kurtosis_VF_V_scatter_blackandwhite_labels_280115.tiff", units="mm", width=xxx, height=xxx, res=300)
plot(varFreq_dVar,log10(kurt_by_mod),bg=c("black","black","black","black","black","white","white","white"),cex=2,pch=21,xlab="variance weighted mean frequency",ylab="log_10(mean kurtosis of trait distribution)")
text(varFreq_dVar[1],log10(kurt_by_mod[1]),"model a: strong selection, tau=1",cex=0.7,pos=4)

text(varFreq_dVar[2],log10(kurt_by_mod[2]),"model b: weak selection, tau=1",cex=0.7,pos=2)

text(varFreq_dVar[3],log10(kurt_by_mod[3]),"model c: strong selection, tau=0",cex=0.7,pos=2)

text(varFreq_dVar[4],log10(kurt_by_mod[4]),"model d: strong selection, tau=0.5",cex=0.7,pos=2)

text(varFreq_dVar[5],log10(kurt_by_mod[5]),"model e: strong selection, tau=2",cex=0.7,pos=4)

text(varFreq_dVar[6],log10(kurt_by_mod[6]),"model f: strong selection, A=B=10000",cex=0.7,pos=4)

text(varFreq_dVar[7],log10(kurt_by_mod[7]),"model g: strong selection, A=4000, B=40",cex=0.7,pos=2)

text(varFreq_dVar[8],log10(kurt_by_mod[8]),"model h: strong selection, A=B=4000",cex=0.7,pos=4)

dev.off()

