setClass("params", slots=list(project_dir = "character", dataset_dir="character", dataset_name="character",
                              seed="numeric", p="numeric", eta="numeric", cancer_pattern_num="numeric",
                              healthy_pattern_num="numeric", convergence_threshold="numeric",
                              maximum_iterations="numeric", prior="character", class_num="numeric",
                              train_sample_num="vector", test_sample_num="vector",
                              train_gt_label="vector", test_gt_label="vector",
                              oversample="logical", evaluate_deconvolution="logical",
                              ratio="numeric", knn="numeric"))

param <- new("params",
             project_dir="~/Documents/Dogs/SRFD-Bayes",
             dataset_dir='simulation_dataset', # choose the dataset type: 'simulation_dataset' and 'real_dataset'
             dataset_name='CNV10', # choose: simulation_dataset = CNV10, CNV30, CNV50; real_dataset = validation_real_data (GSE122126, GSE108462, GSE129373), Xu_data, and Chen_data
             seed=666,
             p=0.5,
             eta=1000,
             cancer_pattern_num=2,
             healthy_pattern_num=7,
             convergence_threshold=1e-2,
             maximum_iterations=1000,
             prior="svm",
             class_num=6,
             oversample=FALSE,
             evaluate_deconvolution=FALSE,
             ratio=1,
             knn=10)

set.seed(param@seed)

if (param@dataset_dir == "simulation_dataset") {
    save_path <- paste(param@project_dir, "/results/", param@dataset_dir, sep="")
} else {
    save_path <- paste(param@project_dir, "/results/", param@dataset_dir, "/", param@dataset_name, sep="")
}
dir.create(save_path, recursive = T)

## load training and  test data
file_dir <- paste(param@project_dir, "/data/", param@dataset_dir, sep="")

if (param@dataset_dir == "simulation_dataset") {

    cat("Dataset: simulation dataset\n")
    cat("CNV event probability <- ", param@dataset_name, "\n")
    train.data <- readMat(paste(file_dir, "/train_data.mat", sep=""))$train.data
    test.data <- readMat(paste(file_dir, "/", param@dataset_name, "/test_data.mat", sep=""))$test.data
    train.theta <- readMat(paste(file_dir, "/train_theta.mat", sep=""))$train.theta
    test.theta <- readMat(paste(file_dir, "/", param@dataset_name, "/test_theta.mat", sep=""))$test.theta    
    param@evaluate_deconvolution <- FALSE
    return(list(train.data, test.data, train.theta, test.theta))
    
} else if (param@dataset_dir == "real_dataset"){
    cat("Dataset: real dataset\n")
    if (param@dataset_name == "validation_real_data"){
        validation.real.data <- readMat(paste(file_dir, "/validation_real_data.mat", sep=""))$validation.real.data
        train.reference <- readMat(paste(file_dir, "/train_reference.mat", sep=""))$train.reference
        return(list(validation.real.data, train.reference))
    } else {
        train.data <- readMat(paste(file_dir, "/", param@dataset_name, "/train_data.mat", sep=""))$train.data
        test.data <- readMat(paste(file_dir, "/", param@dataset_name, "/test_data.mat", sep=""))$test.data
        param@evaluate_deconvolution <- FALSE
        if(param@dataset_name=="Chen_data"){
            param@oversample <- TRUE
        } else {
            param@oversample <- FALSE
        }
        return(list(train.data, test.data, train.theta=NULL, test.theta=NULL))
    }
}

