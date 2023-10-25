#' Get marker genes
#'
#' From a logistic regression model (i.e., the output of trainModel), extract the vector of non-zero coefficients.
#'
#' @param fit The model fit.  Can be either a single glmnet fit object or a named list of them.
#' @param lambda1se Use 1 s.e. Lambda instead of minimum. 
#' @param gns A data.frame to convert the gene given to symbols. The first column must match the rownames of the fit and the second the symbol you want to use.


getMarkers = function(fit,lambda1se=TRUE,gns=NULL){
  if(class(fit)=='cv.glmnet'){
    fit = list(fit=fit)
  }
  out=list()
  for(nom in names(fit)){
    tgt = ifelse(lambda1se,fit[[nom]]$lambda.1se,fit[[nom]]$lambda.min)
    w = which.min(abs(fit[[nom]]$lambda-tgt))
    x = fit[[nom]]$glmnet.fit$beta[,w]
    x = x[x!=0]
    out[[nom]] = data.frame(gene = names(x),
                            coef = x,
                            class = rep(nom,length(x)))
  }
  out = do.call(rbind,out)
  out = out[order(out$class,abs(out$coef)),]
  rownames(out) = NULL
  if(!is.null(gns))
    out$symb = gns[match(out$gene,gns[,1]),2]
  return(out)
}


gm <- readRDS('/icgc/dkfzlsdf/analysis/OE0219_projects/JMMLC/scRNA_Result/LogisticRegression/REF.TrainModel.rds')
GeneIDMap <- read.table('genes.tsv',header = F,stringsAsFactors = F)
genelists = getMarkers(gm,gns=GeneIDMap)
