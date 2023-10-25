#' Code to do the cell matching from Young et al., Science, 2018.
#'
#' Train a regularised logistic regression model on some training single cell data, usually reference normal.  Then use this trained model to infer the similarity of every cell in some target data set to the reference.  A basic plotting function to visualise the results.
#'
#' General notes:
#'
#' The first thing you need to do is load this code by running source('<Path_to_this_file>').
#'
#' The basic data you need for this to be useful are two single cell data-sets (I've assumed 10X, but in theory anything should work).  One of those data sets must be your 'reference' or 'training' data-set that you want to match the second data set (the 'test' or 'target') to.  The training data-set must have an annotation; that is each cell must have a class associated with it.  This doesn't have to be anything special, it could be the cluster number if you want to do something quickliy.
#'
#' The idea is that you first train the logistic regression model using the 'trainModel' function on your training data-set.  This will produce a fit object, or more accurately a list of fit objects, one per class in your annotation.  The main thing to think about when performing this step is how to label each cell into classes.  The way the model fit works, each class is considered one at a time and genes are identified that can be used to tell it apart from all other cells in the data set.  So if you have multiple classes that really represent the same thing (e.g., you've over clustered your data), then the model will be confused as it will see things in the "out" group that really should be in the "in" group.
#'
#' Sometimes the fit will fail for a particular class.  This is usually because there are too few cells in a class to get a meaningful fit.  In this case you can try passing the trainModel nfold=x, where x is some number less than 10 (which is the default).  However, I wouldn't recommend it, as any fit will be of dubious biological value.  When a fit fails for one class, the whole function will still proceed, it will just infer NA values for that class when you apply it later on.
#'
#' Fitting the model can take a loooooong time on even a moderate amount of data (>10,000 cells).  My current approach to this is to intelligently sub-sample the data in such cases, which will happen automatically.  But you could also choose to be very patient or implement something more clever.
#'
#' Having trained the model, you then pass the resulting object and your second single cell data set (the 'test' or 'target' data) to the predictSimilarity function.  This will return logit similarity scores for each cell by default.  If you have some annotation or clustering of the target data set, you can supply it to the predictSimilarity function via the "classes" parameter and it will instead give you the predicted similarity on a per-class level.  This is what you will want most of the time.
#'
#' Finally, there is a simple little plotting function which can aid in visualising the results.  This is called similarityHeatmap.

#############
# Libraries #
#############
library(glmnet)
if(!suppressWarnings(require(ComplexHeatmap)))
  message("ComplexHeatmap not found.  This library is needed for plotting functions.  It's the best heatmap library, so install it anyway :)")
if(!suppressWarnings(require(doMC)))
  message("doMC not found. Parallel processing will not work without this library.")

#############
# Functions #
#############

#' Train the model against a reference
#' 
#' Fit the logistic regression model with regularisation using a One versus Rest approach for the single cell data given in the matrix ref, with cluster labels given.  As the model fit procedure is slow, the number of cells is downsampled if it exceeds some threshold.  To preserve the signal from smaller cell clusters, this downsampling is done in such a way as smaller clusters are not downsampled below a certain threshold.
#'
#' @param refMat A spares matrix, with genes as rows, cells as columns on which the model is to be trained.  Data should be in raw format (i.e., counts).
#' @param classes Vector whose length must match the number of columns in \code{ref}.  Each entry indicates which cluster a particular cell belongs to.
#' @param maxCells Downsample the reference data intelligently if the dataset contains more than this many cells.
#' @param minCells When downsampling, do not let number of cells in a cluster be reduced below this number due to downsampling.
#' @param ... Extra parameters passed to fitting function.
#' @return A model fit, to be used to predict similarities on the target dataset.
trainModel = function(refMat,classes,maxCells=2000,minCells=100, nfolds = 10,nParallel = 20, ...){
  ##############
  # Sub-sample #
  ##############
  nCells = length(classes)
  if(nCells<=maxCells){
    w = seq_along(classes)
  }else{
    cnts = table(classes)
    #Sub-sampling ensures those below minCells are untouched and no group goes below minCells after sampling.  
    ssRate = minCells/cnts
    ssRate[ssRate>1] = 1
    #Really inefficient root finder, but should work reliablyish
    optFun = function(e) abs(sum(cnts*((1-ssRate)*e+ssRate))-maxCells)
    x = optimize(optFun,c(0,1))$minimum
    ssRate = (1-ssRate)*x+ssRate
    w = split(seq_along(classes),classes)
    w = unlist(lapply(names(ssRate),function(e) sample(w[[e]],ceiling(length(w[[e]])*ssRate[e]))))
    w = sort(w)
    message(sprintf("Sub-sampled from %d cells to %d cells.",nCells,length(w)))
  }
  ###################
  # Train the model #
  ###################
  #Prepare the data
  dat = t(refMat[,w])
  dat = dat/Matrix::rowSums(dat)*1e4
  dat = log(dat+1)
  #Train the model
  fit = multinomialFitCV(dat,classes[w],nParallel = nParallel, nfolds = nfolds)
  #####################
  # Return the result #
  #####################
  return(fit)
}

#' Use trained model to predict similarity on target data
#'
#' Given the output of \code{trainModel} and a target single cell data-set, calculate the similarity score for each cluster.  Optionally can merge estimates per-cluster.  Classes that could not be trained in the reference data are given NA similarity.
#'
#' @param fit Model fit on reference data, from \code{trainModel}
#' @param tgtData Sparse single cell reference data matrix.  Columns are cells, rows are genes.  Gene names must match gene names use when training reference data.  Data should be in raw format (i.e., counts).
#' @param classes If not NULL, merges similarities to give a per-class/cluster similarity score.
#' @param minGeneMatch If the fraction of matching genes is less than this, do not proceed.
#' @param logits If FALSE, convert logits into probabilities.
#' @return A similarity matrix giving similarity scores for cells or clusters (columns) against each population in the reference data set(rows).
predictSimilarity = function(fit,tgtData,classes=NULL,minGeneMatch=0.99,logits=TRUE){
  ###################
  # Reguralise data #
  ###################
  tgtGenes = rownames(fit[[1]]$glmnet.fit$beta)
  if(!all(tgtGenes %in% rownames(tgtData))){
    warning(sprintf("Of %d genes used in training data, only %d found in current data set.",length(tgtGenes),sum(tgtGenes %in% rownames(tgtData))))
    if(sum(tgtGenes %in% rownames(tgtData))/length(tgtGenes) < minGeneMatch){
      stop("Could not match enough genes between training data and tgtData")
    }
    tgtGenes = tgtGenes[tgtGenes %in% rownames(tgtData)]
  }
  #Order them as in reference data
  tgtData = tgtData[tgtGenes,]
  #Process in same way
  dat = t(tgtData)
  dat = dat/Matrix::rowSums(dat)*1e4
  dat = log(dat+1)
  ######################
  # Infer similarities #
  ######################
  preds = list()
  for(mark in names(fit)){
    message(sprintf("Predicting similarities for class %s",mark))
    #Do the prediction manually so we can see what is happening
    #Get the lambda value to use
    m = match(fit[[mark]]$lambda.1se,fit[[mark]]$glmnet.fit$lambda)
    #Is there no good match?  If so, set NAs
    if(is.null(fit[[mark]]) | m==1){
      preds[[mark]] = data.frame(mark = rep(NA,nrow(dat)))
    }else{
      preds[[mark]] = data.frame(mark = Matrix::colSums(t(dat) * fit[[mark]]$glmnet.fit$beta[tgtGenes,m]))
    }
    #Name the column
    colnames(preds[[mark]]) = mark
    #preds[[mark]] = predict(fit[[mark]],newx=dat[,rownames(fit[[mark]]$glmnet.fit$beta)],s='lambda.1se',newoffset=rep(0,nrow(dat)))
  }
  pp = do.call(cbind,preds)
  ################
  # Post-process #
  ################
  #Return cluster level summary
  if(!is.null(classes)){
    #If averaging across cluster, this **MUST** be done on probabilities not logits (logits are unbounded and so vulnerable to skewing the answer)
    pp = (1+exp(-pp))**-1
    pp = apply(pp,2,function(e) sapply(split(e,classes),mean))
    #If logits are then required, convert probability back to logits after the averaging using probabilities
    if(logits)
      pp = log(pp)-log(1-pp)
  }else{
    #For non-averaged things, it's up to you which you want.
    if(!logits)
      pp = (1+exp(-pp))**-1
  }
  return(pp)
}

#' Results heatmap
#'
#' Create a heatmap showing the results.  Works best when \code{classes} is supplied to \code{predictSimilarity).
#'
#' @param sims Similarity scores as calculated by \code{predictSimilarity}.
#' @param ... Extra parameters to pass to ComplexHeatmaps
#' @return Heatmap object as constructed by Heatmap
similarityHeatmap = function(sims,isLogit = FALSE, ...){
  #Convert to matrix
  sims = as.matrix(sims)
  #Set colour scheme differently for logits and probabilities
  #isLogit = any(sims<0) | any(sims>1)
  probCols = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
  logitCols = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
  if(isLogit){
    cols = circlize::colorRamp2(seq(-5,5,length.out=length(logitCols)),logitCols)
  }else{
    cols = circlize::colorRamp2(seq(0,1,length.out=3),c("grey","white","red"))
  }
  #Do we have too many target "entities" to show them individually?
  showTgts = nrow(sims)>50
  hm = Heatmap(sims,
               name=paste0('Predicted\nSimilarity\n',ifelse(isLogit,'(Logit)','(Prob)')),
               col = cols,
               show_row_names = !showTgts,
               cluster_rows = showTgts,
               row_title = 'Target Data',
               show_column_names = TRUE,
               cluster_columns = FALSE,
               column_title = 'Reference classes',
               ...
  )
  return(hm)
}


#' Fit the model
#'
#' Do the OvR fit for every variable.  This just does a simple cross-validation selection of regularisation amount.  Trys to handle common failure modes gracefully.
#'
#' @param x Data matrix
#' @param y Class definitions
#' @param nParallel How many parallel processes to do.  Uses DoMC.
#' @param ... Passed to cv.glmnet
#' @return A list containing a model fit object for each unique class in \code{y} or NULL if no model could be fit.
multinomialFitCV = function(x,y,nParallel=1,regAlpha = 0.99, nfolds = 10,...){
  fits = list()
  message(sprintf("nfolds= %s",nfolds))
  if(nParallel>1)
    registerDoMC(cores=nParallel)
  #Do them in order of size
  marks = names(sort(table(as.character(y))))
  for(mark in marks){
    message(sprintf("Fitting model for variable %s",mark))
    fac = factor(y==mark)
    #The two main modes of failure are too few positives and errors constructing lambda.  These should be handled semi-gracefully
    fits[[mark]] = tryCatch(
      cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',nfolds = nfolds,intercept=FALSE,alpha=0.99,type.measure='class',parallel=nParallel>1,...),
      error = function(e) {return(NULL)}
    )
    #A bad fit can either fail entirely, or select a fit with no coefficients
    m = match(fits[[mark]]$lambda.1se,fits[[mark]]$glmnet.fit$lambda)
    #Deviance allows more fine-grained selections to be made when class error is too coarse.  Can be particularly useful when above model defaults to all-zero coefficients
    if(is.null(fits[[mark]]) || is.null(m) || m==1 ){
      if(is.null(m)){
        warning("Default model set all coefficients to zero, refitting using deviance as CV-measure.")
      }else{
        warning("No fit found, refitting using deviance as CV-measure.")
      }
      fits[[mark]] = tryCatch(
        cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',nfolds = nfolds,intercept=FALSE,alpha=0.99, type.measure='deviance',parallel=nParallel>1,...),
        error = function(e) {
          warning(sprintf("Could not fit model for variable %s",mark))
          return(NULL)
        })
    }
    #If we still have an "m=1" error, just return the model, it can be handled at the user's discression at the prediction stage.
  }
  return(fits)
}

#' Get offset for each population/cluster
#'
#' Calculate appropriate intercept term to make the training insensitive to the observed frequency of different populations.
#'
#' @param y A factor or string with two-levels to calculate the appropriate offset.
#' @return A vector of offsets to use in the model fit.
getPopulationOffset = function(y){
  if(!is.factor(y))
    y=factor(y)
  if(length(levels(y))!=2)
    stop("y must be a two-level factor")
  off = sum(y==levels(y)[2])/length(y)
  off = log(off/(1-off))
  return(rep(off,length(y)))
}
