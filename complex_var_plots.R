#Project: Moran modelling RT2
#Creator: Teri Davies
#Date created: 06/03/13
#Date last updated:11/01/14
#Date last run: 11/01/14
#Script purpose: to manipulate the trait_effects.txt file and re-create the EW(2010) plots - part i: individual replicate analysis

args <- commandArgs(TRUE)

trait_filename<-args[1]
filerows<-as.numeric(args[2])
filecols<-as.numeric(args[3])
mut_filename1<-args[4]
filerows1<-as.numeric(args[5])
filecols1<-as.numeric(args[6])
pop<-as.numeric(args[7])
run<-args[8]

#read in the trait effect data
data<-matrix(scan(file=trait_filename, what="raw",n=filerows*filecols,sep="\t"),filerows,filecols,byrow=TRUE)
#make the data numeric 
mode(data)="numeric"
#make a matrix of individual by trait value
trait_val=cbind(data[,1],data[,3])
#make a matrix of individual by mutational effects
#remove the trait column
mut_val=data[,-3]
#remove the fitness column
mut_val=mut_val[,-2]
#make the first column the row names
rownames(mut_val)=mut_val[,1]
mut_val=mut_val[,-1]
message("trait effect data processed")

#read in the mutation freq data 
data1<-matrix(scan(file=mut_filename1, what="raw",n=filerows1*filecols1,sep="\t"),filerows1,filecols1,byrow=TRUE)
#make a vector of the mutation frequencies. The mutations should be in order of age and only those which are not extinct or present in everyone are included
mut_freq=as.numeric(data1[,2])
#transpose the frequencies to make a row vector
mut_freq=t(mut_freq)
#add this mutation frequency row to the mutational effects matrix
mut_val=rbind(mut_freq,mut_val)
message("frequency data processed")

#calculating the trait variance explained per allele 
var_sum=var(trait_val[,2])
message("the sample variance in the trait is ", var_sum)
totvarfile=paste("total_var_",run,".txt",sep="")
write(var_sum,totvarfile)


mut_val_sub=mut_val[2:nrow(mut_val),]

#assign each column a "covariance value" calculated as the variance of that column, plus sum(cov(single mut,all other muts))

#begin by adding a covariance row to the mut_val matrix
covar_vec=c(1:ncol(mut_val))
for (j in 1:length(covar_vec)){
	covar_vec[j]=0
	}
mut_val=rbind(as.numeric(covar_vec),mut_val)

#calculate the covariance matrix of mut_val. Each (i,j) entry is the covariance between the ith and jth columns
#the sum of each column of the covariance matrix is the contribution to the variance in the trait for that mutation
cov_mut_val=cov(mut_val_sub)
message("covariance matrix calculated")

for (i in 1:ncol(mut_val)){
covs=cov_mut_val[,i]
covs=sort(covs)
mut_val[1,i]=covs[1]
for (j in 2:length(covs)){
mut_val[1,i]=mut_val[1,i]+covs[j]
}
message("mutation ",i," complete")
}

#check that the mutation specific variances sum to the total variance
checkit=mut_val[1,]
checkit=sort(checkit)
checksum=checkit[1]
for (i in 2:length(checkit)){
checksum=checksum+checkit[i]
}

message("the allele variances sum to ", checksum)


#pool the variances of alleles with the same frequency
var_plot=matrix(NA,nrow=(pop-1),2) 

#assign the first column of the plot matrix to be these frequencies (freq=1...(pop-1))
var_plot[,1]=c(1:(pop-1))
message("frequencies assigned to column 1")

new_mut_val=mut_val[1:2,]
#sort the covariances and frequencies by covariance so that we begin with the smallest covs for precision reasons
new_mut_val=new_mut_val[,sort.list(new_mut_val[1,])]

#now populate the covariance column of var_plot by pooling variances for alleles with the same frequency
for (j in 1:length(new_mut_val[1,])){ #for each mutation
for (i in 1:length(var_plot[,1])){ # for each potential freq
if (var_plot[i,1]==new_mut_val[2,j]){ # if this freq=the freq for our current mutation
if (is.na(var_plot[i,2])==TRUE){
var_plot[i,2]=new_mut_val[1,j]
}
else {
var_plot[i,2]=var_plot[i,2]+new_mut_val[1,j] # update the variance for this freq
}
message("variance for mutation ",j," contributed to matrix")
break # break out of the for loop and move on to the next mutation
}
}
}


#write the allele frequencies plus covariance sums to file for pooling across replicate runs
outfile=paste("var_plot_",run,".txt",sep="")
write.table(var_plot,outfile,sep="\t",eol="\n",quote=FALSE,row.names=FALSE, col.names=FALSE)

#check that the trait value for each individual equals the sum of the mutational trait effects

trait_check=matrix(NA,nrow=filerows,ncol=3) 

for (i in 1:filerows){
trait_check[i,1]=data[i,3] # the trait value outputted by the c code
traits=sort(mut_val_sub[i,])
trait_check[i,2]=traits[1]
for (j in 2:(filecols-3)){
trait_check[i,2]=trait_check[i,2]+traits[j] # the calculated sum of trait effects after ordering the effects
}
}

trait_check[,3]=trait_check[,1]-trait_check[,2]

outfile=paste("trait_check_",run,".txt",sep="")
write.table(trait_check,outfile,sep="\t",eol="\n",quote=FALSE,row.names=FALSE, col.names=FALSE)




