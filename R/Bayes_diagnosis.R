#x = train_data
#y = param@train_gt_label
#newx = test_data
#prior = param@prior
#sf = list(train = train_sf, test = test_sf)
#method = "SVM"
#family = "binomial"
#train_V_without_label = train_sf
#test_V_without_label = test_sf

# Presuming this is the same as betafit
estBetaParams <- function(x) {
  mu <- mean(x)
  var <- var(x)
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

Bayes_diagnosis <- function(x, y, newx, method="SVM", train_V_without_label, test_V_without_label, param){
#  library(brms) # install from github (https://github.com/paul-buerkner/brms)
  if (method == "rf") { # Random forest
    library(caret)
    fitControl <- trainControl(method="oob", number=100, classProbs=T, savePred=T)
    df.x <- data.frame(t(x))
    df.y <- data.frame(y)
    df.y$y <- factor(df.y$y, levels = unique(df.y$y), labels=paste0("C", unique(df.y$y)))
    rf_mod <- train(x=df.x, y=df.y$y,  method="rf", norm.votes=T, type="Classification", ntree=100, trControl=fitControl)
    posterior <- predict(rf_mod, newdata=data.frame(t(newx)), type="prob")
    pre_results <- as.character(predict(rf_mod, newdata=data.frame(t(newx))))
    pre_labels <- as.numeric(gsub("C", "", pre_results))
    rm(df.x, df.y, rf_mod)
    
  } else if (method == "svm") { # support vector machine 
    library(e1071)
    mod.x <- cbind(as.factor(y), t(x))
    colnames(mod.x) <- c("y", seq(1,ncol(mod.x)-1))
    svm.model <- svm(y~., mod.x, type = "C-classification", kernel = "linear", scale = TRUE, cross = 10, probability = TRUE)
    new.x <- t(newx)
    colnames(new.x) <- seq(1,ncol(new.x))
    svm.predict <- predict(svm.model, new.x, type="response", probability = TRUE)
    pre_labels <- as.numeric(svm.predict)
    posterior <- attr(svm.predict, "probabilities")
    
  } else if (method == "knn") { # k-nearest neighbour
    library(caret)
    #preProcess <- c("center", "scale")
    fitControl <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=T, savePred=T)
    df.x <- data.frame(t(x))
    df.y <- data.frame(y)
    df.y$y <- factor(df.y$y, levels = unique(df.y$y), labels=paste0("C", unique(df.y$y)))
    #knn_mod <- train(x=df.x, y=df.y$y,  method="knn", metric='Accuracy', preProcess = preProcess, trControl=fitControl)
    knn_mod <- train(x=df.x, y=df.y$y,  method="knn", metric='Accuracy', trControl=fitControl)
    posterior <- predict(knn_mod, newdata=data.frame(t(newx)), type="prob")
    pre_results <- as.character(predict(knn_mod, newdata=data.frame(t(newx))))
    pre_labels <- as.numeric(gsub("C", "", pre_results))
    rm(df.x, df.y, knn_mod)
    
  } else if (method == "nnet") { # neural network
    library(caret)
    fitControl <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=T, savePred=T)
    df.x <- data.frame(t(x))
    df.y <- data.frame(y)
    df.y$y <- factor(df.y$y, levels = unique(df.y$y), labels=paste0("C", unique(df.y$y)))
    nnet_mod <- train(x=df.x, y=df.y$y,  method="nnet", linout=TRUE, trControl=fitControl, trace=FALSE)
    posterior <- predict(nnet_mod, newdata=data.frame(t(newx)), type="prob")
    pre_results <- as.character(predict(nnet_mod, newdata=data.frame(t(newx))))
    pre_labels <- as.numeric(gsub("C", "", pre_results))
    rm(df.x, df.y, nnet_mod)
    
  } else if (method == "nb") { # naive bayes (try nb discrete - nbDiscrete)
    library(caret)
    library(klaR) # for nb
    #library(bnclassify) # for nbDiscrete
    fitControl <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=T, savePred=T)
    df.x <- data.frame(t(x))
    df.y <- data.frame(y)
    df.y$y <- factor(df.y$y, levels = unique(df.y$y), labels=paste0("C", unique(df.y$y)))
    nb_mod <- train(x=df.x, y=df.y$y,  method="nb", trControl=fitControl)
    posterior <- predict(nb_mod, newdata=data.frame(t(newx)), type="prob")
    pre_results <- as.character(predict(nb_mod, newdata=data.frame(t(newx))))
    pre_labels <- as.numeric(gsub("C", "", pre_results))
    rm(df.x, df.y, nb_mod)

  } else if (method == "mlp") { # multi-layer perceptron
    library(caret)
    library(RSNSS)
    fitControl <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=T, savePred=T)
    df.x <- data.frame(t(x))
    df.y <- data.frame(y)
    df.y$y <- factor(df.y$y, levels = unique(df.y$y), labels=paste0("C", unique(df.y$y)))
    nb_mod <- train(x=df.x, y=df.y$y,  method="mlp", trControl=fitControl)
    posterior <- predict(nb_mod, newdata=data.frame(t(newx)), type="prob")
    pre_results <- as.character(predict(nb_mod, newdata=data.frame(t(newx))))
    pre_labels <- as.numeric(gsub("C", "", pre_results))
    rm(df.x, df.y, nb_mod)
    
  }
  
  # Bayes
  class_type <- unique(y)
  theta_dist <- list()
  for (i in 1:length(class_type)){
    train_burden <- train_V_without_label[(param@healthy_pattern_num+1):nrow(train_V_without_label), y == class_type[i]]
    for (j in 1:nrow(train_burden)){
      train_theta_one <- train_burden[j,]
      beta_ab <- estBetaParams(train_theta_one)
      if(j==1) {
        theta_dist[[i]] <- data.frame(beta_ab)
      } else {
        theta_dist[[i]] <- rbind(theta_dist[[i]], data.frame(beta_ab))  
      }
    }
  }
  bayes_prob <- t(posterior)
  test_burden <- test_V_without_label[(param@healthy_pattern_num + 1):nrow(test_V_without_label),]
  for (i in 1:ncol(test_V_without_label)) {
    for (j in 1:length(class_type)) {
      for (k in 1:nrow(theta_dist[[j]])) {
        lik <- dbeta(test_burden[k,i], theta_dist[[j]][k,]$alpha, theta_dist[[j]][k,]$beta)
        if(lik>0 & is.infinite(lik)==F) bayes_prob[j,i] <- bayes_prob[j,i] * lik # had to add test here because of 0 or Inf values
      }
    }
    bayes_prob[,i] <- bayes_prob[,i]/sum(bayes_prob[,i])
  }
  bayes_labels <- apply(bayes_prob, 2, which.max)
  return(list(posterior=posterior, pre_labels=pre_labels, bayes_prob=bayes_prob, bayes_labels=bayes_labels))
}
