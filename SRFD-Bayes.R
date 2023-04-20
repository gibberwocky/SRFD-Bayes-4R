library(matlab)
library(R.matlab)

## initialization
source("~/Documents/Dogs/SRFD-Bayes/R/Initialization.R")
source("~/Documents/Dogs/SRFD-Bayes/R/deconvolution_with_reference.R")
source("~/Documents/Dogs/SRFD-Bayes/R/SRFD.R")
source("~/Documents/Dogs/SRFD-Bayes/R/SRFD_Bayes.R")
source("~/Documents/Dogs/SRFD-Bayes/R/Bayes_diagnosis.R")
source("~/Documents/Dogs/SRFD-Bayes/R/evaluate_deconvolution.R")
source("~/Documents/Dogs/SRFD-Bayes/R/evaluate_diagnosis.R")
source("~/Documents/Dogs/SRFD-Bayes/R/borderline_smote.R")
source("~/Documents/Dogs/SRFD-Bayes/R/mix_WH.R")

if(param@dataset_dir == "real_dataset" & param@dataset_name == "validation_real_data"){
  validation_tf <-  deconvolution_with_reference(validation.real.data, train.reference, rep(0,ncol(validation.real.data)))
  save(validation_tf, file=paste(save_path, "/validation_tf.RData", sep=""))
} else {
  ## SRFD-Bayes algorithm
  #-------------------------------------------------------------------#
  #----Input:
  #------training data matrix, test data matrix, parameters
  #----Output:
  #------diagnostic results:
  #------The first row represents the predicted class clabels
  #------The second row represents the predicted tumor fraction
  #[train_reference, SRFD_Bayes_results] <- SRFD_Bayes(train.data, test.data)
  SRFD_Bayes_results <- SRFD_Bayes(train.data, test.data, method="svm")
  #-------------------------------------------------------------------#
  
  ## save results
  #save(strcat(save_path,'/SRFD_Bayes_results.mat'), 'SRFD_Bayes_results')
  save(SRFD_Bayes_results, file=paste(save_path, "/SFRD_Bayes_results.RData", sep=""))
  
  if(param@dataset_dir=="simulation_dataset") {
    #save(strcat('../data/real_dataset/train_reference.mat'), 'train_reference') 
    ## evaluate deconvolution (simulation dataset only)
    if (evaluate_deconvolution) {
      cat("\nEvaluating the deconvolution performance\n")
      evaluate_deconvolution(SRFD_Bayes_results$tumor_fraction, test.theta)
    }
  }
  
  ## evaluate diagnosis
  if (param@class_num == 2){ # Not yet implemented
      #[sorted_prob, sorted_index] <- sort(SRFD_Bayes_results(1,:))
      #sorted_label <- param.test_gt_label(sorted_index)-1
      #[FPR,TPR] <- cal_roc(sorted_prob, sorted_label)
      #auc <- abs(trapz(FPR,TPR))
      #disp(['AUC <- ' num2str(auc)])
  } else {
      cat("\n>>>>>>>>Evaluating the diagnostic performance...\n")
      average_ACC <- evaluate_diagnosis(SRFD_Bayes_results)
      cat("\naverage ACC: ", average_ACC)
  }
}

# Chen data
# svm     average ACC:  0.5046778
# rf      average ACC:  0.5340333
# knn     average ACC:  0.4045839 (slow)
# nnet    average ACC:  0.3543767
# nb      average ACC:  0.2177545 (very slow)
# mlp     average ACC:  0.4473385

# >>>>>>>>Evaluating the diagnostic performance... (Chen, svm)
#confusion_matrix: 
#  [,1] [,2] [,3] [,4] [,5] [,6]
#[1,]  203    0    0    1    1    2
#[2,]    2    0    0    0    0    0
#[3,]    0    0    0    0    2    1
#[4,]    6    0    0   22    1    0
#[5,]    2    0    0    0    9    0
#[6,]   13    0    0    0    5   16
#average ACC:  0.5046778

# >>>>>>>>Evaluating the diagnostic performance... (Xu, svm)
#confusion_matrix: 
#  [,1] [,2]
#[1,]  368   49
#[2,]   58  288
#average ACC:  0.857432
 
# >>>>>>>>Evaluating the diagnostic performance... (CNV30, svm)
#confusion_matrix: 
#  [,1] [,2] [,3] [,4] [,5] [,6]
#[1,]  328    3   20   27   18    4
#[2,]   14  368    0   13    1    4
#[3,]   18    2  374    5    1    0
#[4,]   14    7    2  372    4    1
#[5,]   22    0    0    1  377    0
#[6,]   10    1    0    6    0  383
#average ACC:  0.9175

