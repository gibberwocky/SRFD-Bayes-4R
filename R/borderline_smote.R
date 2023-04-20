#train_R <- x$R
#train_labels <- param@train_gt_label
#post_num <- 40
#data_type <- "fraction"

compute_tumor_dis <- function(single_h, train_R, param, data_type) {
  if (data_type == "fraction") {
    tumor_h <- pracma::kron(diag(param@class_num-1), matrix(1, nrow=1, ncol=param@cancer_pattern_num)) %*% single_h[(param@healthy_pattern_num+1):length(single_h),]
    tumor_train <- kronecker(diag(param@class_num-1), matrix(1, nrow=1, ncol=param@cancer_pattern_num)) %*% train_R[(param@healthy_pattern_num+1):nrow(train_R),]
    dis_out <- colSums((tumor_h %*% t(matrix(1, ncol(tumor_train), 1)) - tumor_train)^2)
  } else if (data_type == "methylation") {
    dis_out <- colSums((single_h %*% t(single_h)) - 2 * single_h %*% t(train_R) + (train_R %*% t(train_R)))
  }
  return(dis_out)
}


borderline_smote <- function(train_R, train_labels, post_num,data_type, param) {

  class_index <- seq(1:param@class_num)
  synthetic_data <- NULL
  synthetic_labels <- NULL
  train_data_without_noise <- NULL
  train_labels_without_noise <- NULL
  danger_data <- NULL
  danger_labels <- NULL
  
  for(i in 1:length(class_index)){
    class_data <- train_R[,train_labels == class_index[i]]
    randindex <- sample(ncol(class_data))
    class_data <- class_data[,randindex]
    nearest_in_num <- rep(0, ncol(class_data))
    for (j in 1:ncol(class_data)) {
      single_h <- as.matrix(class_data[,j])
      dis_out <- compute_tumor_dis(single_h, train_R, param, data_type)
      dis_out[dis_out==0] <- Inf
      sorted_index <- as.numeric(sort.list(dis_out))
      nearest_labels <- train_labels[sorted_index[1:param@knn]]
      nearest_in_num[j] <- sum(nearest_labels == class_index[i], na.rm=T)
    }
    
    nearest_in_num[nearest_in_num == 0] <- Inf  # noise data
    names(nearest_in_num) <- seq(1,length(nearest_in_num))
    sorted_num <- as.numeric(sort(nearest_in_num))
    sorted_index <- as.numeric(names(sort(nearest_in_num)))
    sorted_index <- which(sorted_num != Inf)
    
    if (length(sorted_index) == 0) {
      next
    }
    
    danger_data <- cbind(danger_data, class_data[,sorted_index[1:min(param@knn, length(sorted_index))]])
    danger_labels <- c(danger_labels, rep(class_index[i], 1, min(param@knn,length(sorted_index))))
    train_data_without_noise <- cbind(train_data_without_noise, class_data[,sorted_index])
    train_labels_without_noise <- c(train_labels_without_noise, i*ones(1,length(sorted_index)))

  }
  
  extend_index <- unique(danger_labels)
  if(length(extend_index) < param@class_num) {
    test <- 1
  }
  
  extend_num <- NULL
  for(i in seq_along(extend_index)) {
    danger_num <- sum(danger_labels == extend_index[i])
    train_num <- sum(train_labels_without_noise == extend_index[i])
    temp <- (max(c(post_num, danger_num)) - train_num) / danger_num
    extend_num <- c(extend_num, temp)
  }
  extend_num[extend_num < 0] <- 0
  synthetic_data <- NULL
  synthetic_labels <- NULL
  
  for(i in length(extend_index)) {
    class_data <- danger_data[,danger_labels == extend_index[i]]
    class_in_data <- train_data_without_noise[,train_labels_without_noise == extend_index[i]]
    for(j in 1:ncol(class_data)) {
      single_h <- as.matrix(class_data[,j])
      dis_in <- compute_tumor_dis(single_h, class_in_data, param, data_type)
      if(data_type == "methylation") {  # avoid same data
        dis_in[dis_in==0] = Inf
        sorted_dis <- sort(dis_in)
        sorted_index <- sort.list(dis_in)
        sorted_index <- sorted_index[sorted_dis != Inf]
      } else if (data_type == "fraction") {
        sorted_dis <- sort(dis_in)
        sorted_index <- sort.list(dis_in)
      }
      max_sample <- min(param@knn, ncol(sorted_index))
      nearest_H <- class_in_data[,sorted_index[1:max_sample]]  #compute the nearest neighbors
      
      if(extend_num[i] == 0) {
        next
      } else if(extend_num[i] < 1) {
        if(runif(1) > extend_num[i]) {
          next
        } else {
          extend_num_i <- 1
        }
      } else {
        extend_num_i <- round(extend_num[i])
      }
      
      for (k in 1:extend_num_i) {
        rand_neighbor <- nearest_H[,sample(max_sample,1)]  #random neighbor
        mixed_coff <- runif(1) #random mixed coefficient
        synthetic_data <- cbind(synthetic_data, single_h + mixed_coff*(rand_neighbor-single_h))
        synthetic_labels <- c(synthetic_labels, extend_index[i])
      }
      
    }
  }
  return(list(synthetic_data=synthetic_data,synthetic_labels=synthetic_labels))
}
  

