###########################################################
# RF-Score_pred.r (all rights reserved)					  #
# 	Author: Dr Pedro J. Ballester                         #
#	Usage:  Read Appendix_A1.doc  		                  #
###########################################################

library(randomForest)

##################################################
# READING TRAINING (TRN) AND TEST (TST) DATASETS #
##################################################

rawtrndata 	= read.csv("PDBbind_refined07-core07.csv", na.strings=c(".", "NA", "", "?"))
rawtstdata 	= read.csv("PDBbind_core07.csv", 		 na.strings=c(".", "NA", "", "?"))
itarget	= 1						# column of target in table
nfeats 	= ncol(rawtrndata) - 1    		# with target (1), but without pdb code (last column)
ntrndata 	= nrow(rawtrndata)                  # number of pdb complexes for training
ntstdata 	= nrow(rawtstdata) 			# number of pdb complexes for testing
seed		= 1			

itrain	= seq(1,ntrndata,1)
nsample	= ntrndata			
set.seed(seed)		
itrain	= sample(itrain)[1:nsample]		# shuffle selected complexes
trainTarget = rawtrndata[ itrain, itarget]
trainFeats  = rawtrndata[ itrain, 2:nfeats]
rownames(trainFeats)[1:nrow(trainFeats)] = as.vector(rawtrndata[ itrain, nfeats+1])

itest		= seq(1,ntstdata,1)
testTarget  = rawtstdata[ itest, itarget] 
testFeats   = rawtstdata[ itest, 2:nfeats]
rownames(testFeats)[1:nrow(testFeats)]   = as.vector(rawtstdata[ itest, nfeats+1])
	
##################################################
# DATA PRE-PROCESSING					 #
##################################################

ifeats 	= seq(1,ncol(trainFeats),1)    	
selfeats 	= c(1:3,6, 10:12, 15, 19:21, 24, 28:30, 33, 37:39, 42, 46:48, 51, 55:57, 60, 64:66, 69, 73:75, 78) 
nselfeats 	= ncol(trainFeats[,selfeats])       
trainFeats	= trainFeats[,selfeats]	
testFeats	= testFeats[,selfeats]	

##################################################
# SELECTING RF WITH BEST INTERNAL VALIDATION     #
# (RF-SCORE)						 #
##################################################

set.seed(seed)
rmse_OOB_best = 1e10		#dummy high value
for (mtry in 2:nselfeats) 
{
   	RF_mtry = randomForest(trainTarget ~ ., data=trainFeats, 
					   ntree=500, mtry=mtry, na.action=na.omit)
	rmse_OOB = sqrt(mean(( RF_mtry$y - RF_mtry$predicted)^2)) 
      if (rmse_OOB < rmse_OOB_best)
	{ 
		mbest 		= mtry
		rmse_OOB_best 	= rmse_OOB
   		print(paste("mbest=",mbest, "rmse_OOB=",round(rmse_OOB,3)))
	}
	print(paste("mtry=",mtry))
}

RF_Score 	= randomForest(trainTarget ~ ., 
                          data=trainFeats, 
				  ntree=500, mtry=mbest, importance=TRUE, na.action=na.omit) #includes calculation of variable importance

##################################################
# VARIABLE IMPORTANCE BY RF-SCORE			 #
##################################################

RF_Score$importance = importance(RF_Score, type=1) #only %incMSE
windows(height = 15, width = 12); varImpPlot(RF_Score, n.var=nrow(RF_Score$importance), main="Variable Importance")  #default is top 30 features

##################################################
# RF-SCORE FIT OF TRAINING DATASET			 #
##################################################

trainPred 	= predict(RF_Score, newdata = trainFeats) 
rmse 		= format(sqrt(mean(( trainTarget - trainPred)^2)), digits=3)
sdev 		= format(sd( trainTarget - trainPred ), digits=3)
fitpoints 	= na.omit(cbind(trainTarget, trainPred))
fitcorr  	= format(round(cor(fitpoints[,1], fitpoints[,2]), digits=3))      #Pearson correlation coefficient 
sprcorr   	= format(cor(rank(fitpoints[,1]), rank(fitpoints[,2])), digits=3) #Spearman rank correlation 

windows();  plot(fitpoints[,1], fitpoints[,2], asp=1, xlab="Measured binding affinity (PDBbind DB)", ylab="Predicted binding affinity (RF-Score)")
prline 	= lm(fitpoints[,2] ~ fitpoints[,1])						#Linear fit between predicted and measured binding affinity
abline(prline) 
title(main=paste("R=",fitcorr, "on training set (",ntrndata," complexes)"))
grid()  

##################################################
# RF-SCORE PREDICTION OF TEST DATASET	       #
##################################################

testPred 	= predict(RF_Score, newdata = testFeats) 
rmse 		= format(sqrt(mean(( testTarget - testPred)^2)), digits=3) 
sdev 		= format(sd( testTarget - testPred ), digits=3)
fitpoints 	= na.omit(cbind(testTarget, testPred))
fitcorr  	= format(round(cor(fitpoints[,1], fitpoints[,2]), digits=3))   	#Pearson correlation coefficient (R)
sprcorr     = format(cor(rank(fitpoints[,1]), rank(fitpoints[,2])), digits=3) #Spearman rank correlation 

windows();  plot(fitpoints[,1], fitpoints[,2], asp=1, xlab="Measured binding affinity (PDBbind DB)", ylab="Predicted binding affinity (RF-Score)")
prline 	= lm(fitpoints[,2] ~ fitpoints[,1])						#Linear fit between predicted and measured binding affinity
abline(prline) 
title(main=paste("R=",fitcorr, "on independent test set (",ntstdata," complexes)"))
grid()  

##################################################
# WRITING OUT RF-SCORE PREDICTION OF TEST DATASET#
##################################################

write.csv(cbind(testPred, rawtstdata[seq(1,ntstdata,1), nfeats+1]), row.names=rawtstdata[seq(1,ntstdata,1), nfeats+1], file="RF-Score_pred.csv")
