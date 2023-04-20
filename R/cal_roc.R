cal_roc <- function(sorted_tf, sorted_label){

FPR <- zeros(1,size(sorted_label,2))
TPR <- zeros(1,size(sorted_label,2))
for (i in 1:size(sorted_tf,2)){
    pre_label <- [zeros(1,i),ones(1,size(sorted_tf,2)-i)]
    tp <- 0
    fp <- 0
    fn <- 0
    tn <- 0

    for (j in 1:size(pre_label,2)){
        if (sorted_label(j)==1 && pre_label(j)==1){
            tp <- tp+1
        } else if (sorted_label(j)==1 && pre_label(j)==0){
            fn <- fn+1
        } else if (sorted_label(j)==0 && pre_label(j)==1){
            fp <- fp+1
        } else if (sorted_label(j)==0 && pre_label(j)==0){
            tn <- tn+1
        }
    }

    sens <- tp/(tp+fn)
    spec <- tn/(tn+fp)
    TPR(i) <- sens
    FPR(i) <- 1-spec

}
}

