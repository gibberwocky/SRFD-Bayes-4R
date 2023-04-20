#test_tf <- SRFD_Bayes_results$SRFD_Bayes_prediction
#test_theta <- test.theta

evaluate_deconvolution <- function(test_tf,test_theta, param){
  ## healthy evaluation
  cat("\n>>>>>>>>Healthy evaluation<<<<<<<<<<\n")
  hRMSE <- healthy_eval(test_tf,test_theta)
  cat("\nRMSE <- ", hRMSE, "\n")
  ## tumor fraction evaluation
  cat("\n>>>>>>>>Tumor fraction evaluation<<<<<<<<<<\n")
  te <- tumor_eval(test_tf,test_theta)
  cat("\nRMSE <- ", te$RMSE)
  cat("\nPCC <- " , te$PCC)
}

healthy_eval <- function(test_tf, test_theta){
  healthy_theta_pre <- test_tf[1:param@test_sample_num[1]]
  healthy_theta_gt <- 1 - test_theta[1:param@test_sample_num[1]]
  RMSE <- sqrt(sum((healthy_theta_pre - healthy_theta_gt)^2)/length(healthy_theta_pre))
  return(RMSE)
}

tumor_eval <- function(test_tf, test_theta){
  predict_theta <- test_tf[(param@test_sample_num[1]+1):length(test_tf)]
  theta_gt <- test_theta[(param@test_sample_num[1]+1):length(test_theta)]
  RMSE <- sqrt(sum((predict_theta - theta_gt)^2)/length(predict_theta))
  PCC <- cor(predict_theta, theta_gt)
  return(list(RMSE=RMSE,PCC=PCC))
}

