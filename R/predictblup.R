#Function: ghap.predictblup
#License: GPLv3 or later
#Modification date: 14 May 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Predict blup solutions based on reference values

ghap.predictblup <- function(
  refblup,
  vcp,
  covmat,
  errormat = NULL,
  include.ref = TRUE,
  inbreeding = FALSE,
  tol = 1e-10,
){
  
  # Sanity check for input objects ---------------------------------------------
  if(is.null(names(refblup))){
    stop("\nArgument 'refblup' must be a named vector.\n")
  }
  if(identical(rownames(covmat),colnames(covmat)) == FALSE){
    stop("\nRow and column names in 'covmat' are not identical.\n")
  }
  if(length(which(names(refblup) %in% colnames(covmat))) != length(refblup)){
    stop("\nSome 'refblup' names are missing in 'covmat'.\n")
  }
  if(is.null(errormat) == FALSE){
    if(identical(rownames(errormat),colnames(errormat)) == FALSE){
      stop("\nRow and column names in 'errormat' are not identical.\n")
    }
    if(identical(sort(rownames(errormat)),sort(names(refblup))) == FALSE){
      stop("\nNames in 'errormat' must match the names in 'refblup'.\n")
    }
  }

  # Map test and reference -----------------------------------------------------
  ref <- names(refblup)
  test <- colnames(covmat)
  if(include.ref == FALSE){
    test <- test[which(test %in% ref == FALSE)]
  }
  if(inbreeding == FALSE){
    diag(covmat) <- 1    
  }
  
  # Invert reference matrix ----------------------------------------------------
  Krr.i <- try(solve(covmat[ref,ref]), silent = TRUE)
  if(inherits(Krr.i, "try-error")){
    Krr.i <- try(solve(covmat[ref,ref] + Diagonal(length(ref))*tol), silent = TRUE)
    if(inherits(LHSi, "try-error")){
      emsg <- paste0("\nUnable to invert the reference matrix",
                     " even after adding a tolerance of ",
                     tol)
      stop(emsg)
    }
  }
  
  # Solve blup for test individuals --------------------------------------------
  k <- covmat[test,ref]%*%Krr.i
  u <- k%*%refblup
  results <- data.frame(Estimate = u)
  results$`Std. Error` <- NA
  results$`Accuracy` <- NA
  rownames(results) <- test
  
  # Check if errors should be computed -----------------------------------------
  if(is.null(errormat) == FALSE){
    varuhat <- k%*%(vcp*covmat[ref,ref] - errormat)%*%t(k)
    varuhat <- diag(as.matrix(varuhat))
    varu <- vcp*diag(covmat[test,test])
    results$Accuracy <- sqrt(varuhat/varu)
    results$`Std. Error` <- sqrt((1 - results$Accuracy^2)*varu)
  }

  #Return output --------------------------------------------------------------
  return(results)
  
}
