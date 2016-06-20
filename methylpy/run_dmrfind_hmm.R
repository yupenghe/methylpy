library(RHmm)
source("/data/home/mschultz/utilities/my_python_modules/h1_derived/h1_hmm.R")
args = commandArgs(trailingOnly=TRUE)

if(length(args) != 2){
        stop("Usage: Rscript run_dmrfind_hmm.R <chrom number for training> <chrom number for predicting>")
}
print("Loading data for prediction into R")
print(date())
data = read.table(paste("chisquare_results_",args[2],".tsv",sep=""),stringsAsFactors=FALSE)

print("Loading training data into R")
print(date())
train = read.table(paste("chisquare_results_",args[1],".tsv",sep=""),stringsAsFactors=FALSE)
colnames(data) = c("chr","position","strand","class","pvalue","adjusted","test")
colnames(train) = c("chr","position","strand","class","pvalue","adjusted","test")

print("Training HMM")
print(date())
hmm = HMMFit(as.numeric(train$test),nStates=2,dis="DISCRETE",control=list(nInit=100))
#data$adjusted=as.numeric(data$adjusted)
save(hmm,file="hmm.Rd")
print("Running viterbi")
print(date())
results=predict_states(hmm, data,obs_name="test")
write.table(results,file=paste("hmm_results_",args[2],".tsv",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=F)
