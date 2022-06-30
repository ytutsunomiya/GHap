#Function: ghap.simmating
#License: GPLv3 or later
#Modification date: 30 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Simulate individuals from specified matings

ghap.simmating <- function(
  object,
  n.individuals = 1,
  parent1 = NULL,
  parent2 = NULL,
  model = "proportional",
  out.file,
  only.active.markers = TRUE,
  ncores = 1,
  verbose = TRUE
){
  
  # Check if phase is a GHap.phase object-------------------------------------------------------------
  if(inherits(object, "GHap.phase") == FALSE){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Check if inactive markers should be reactivated---------------------------------------------------
  if(only.active.markers == FALSE){
    object$marker.in <- rep(TRUE,times=object$nmarkers)
    object$nmarkers.in <- length(which(object$marker.in))
  }
  
  # Check if out file exist---------------------------------------------------------------------------
  samples.file <- paste(out.file,"samples",sep=".")
  phase.file <- paste(out.file,"phase",sep=".")
  markers.file <- paste(out.file,"markers",sep=".")
  pedigree.file <- paste(out.file,"pedigree",sep=".")
  if(file.exists(samples.file) == TRUE | file.exists(markers.file) == TRUE | file.exists(markers.file) == TRUE){
    stop("Output file already exists!")
  }else{
    rnumb <- runif(n = 1, min = 1, max = 1e+6)
    rnumb <- ceiling(rnumb)
    tmp.file <- paste(tempdir(),"/tmp",rnumb,sep="")
    tmp.samples.file <- paste(tmp.file,"samples",sep=".")
    tmp.pedigree.file <- paste(tmp.file,"pedigree",sep=".")
    tmp.phase.file <- paste(tmp.file,"phase",sep=".")
    tmp.markers.file <- paste(tmp.file,"markers",sep=".")
  }
  
  # Check parent vectors------------------------------------------------------------------------------
  if(is.null(parent1) | is.null(parent2)){
    stop("Missing parents for mating simulations\n")
  }
  if(inherits(parent1, "character")){
    tmpprobs <- rep(1/length(parent1), times = length(parent1))
    tmpnames <- parent1
    parent1 <- tmpprobs
    names(parent1) <- tmpnames
  }
  if(inherits(parent2, "character")){
    tmpprobs <- rep(1/length(parent2),times = length(parent2))
    tmpnames <- parent2
    parent2 <- tmpprobs
    names(parent2) <- tmpnames
  }
  parent1ok <- which(names(parent1) %in% object$id)
  parent2ok <- which(names(parent2) %in% object$id)
  if(length(parent1ok) != length(parent1) | length(parent1ok) != length(parent1)){
    stop("Some of the provided parent names were not found in the phase object\n")
  }
  if(verbose == TRUE){
    cat("\n\nMating simulation started...\n")
    cat("Number of parents in set 1: ", length(parent1), "\n", sep="")
    cat("Number of parents in set 2: ", length(parent2), "\n", sep="")
    cat("Number of progeny: ", n.individuals, "\n\n", sep="")
  }
  
  # Function for individual simulation----------------------------------------------------------------
  indbuild <- function(i){
    crossovers1 <- rpois(n = 1, lambda = chrmean[chr])
    crossovers2 <- rpois(n = 1, lambda = chrmean[chr])
    if(crossovers1 > 0){
      breakpoints1 <- sample(x = 1:m, size = crossovers1,
                             replace = F, prob = probchr)
      breakpoints1 <- sort(breakpoints1)
      parent1idx <- which(colnames(parenthap) == indtbl$parent1[i])
      parent1idx <- sample(parent1idx, size=2, replace=F)
      hap1idx1 <- c(1,breakpoints1+1)
      hap1idx2 <- c(breakpoints1,m)
      hap1idx1 <- hap1idx1[which(1:length(hap1idx1) %% 2 == 0)]
      hap1idx2 <- hap1idx2[which(1:length(hap1idx2) %% 2 == 0)]
      hap1 <- parenthap[,parent1idx[1]]
      for(j in 1:length(hap1idx1)){
        cropseg <- hap1idx1[j]:hap1idx2[j]
        cropseg <- cropseg[which(cropseg < m)]
        hap1[cropseg] <- parenthap[cropseg,parent1idx[2]]
      }
    }else{
      parent1idx <- which(colnames(parenthap) == indtbl$parent1[i])
      parent1idx <- sample(parent1idx, size=2, replace=F)
      hap1 <- parenthap[,parent1idx[1]]
    }
    if(crossovers2 > 0){
      breakpoints2 <- sample(x = 1:m, size = crossovers2,
                             replace = F, prob = probchr)
      breakpoints2 <- sort(breakpoints2)
      parent2idx <- which(colnames(parenthap) == indtbl$parent2[i])
      parent2idx <- sample(parent2idx, size=2, replace=F)
      hap2idx1 <- c(1,breakpoints2+1)
      hap2idx2 <- c(breakpoints2,m)
      hap2idx1 <- hap2idx1[which(1:length(hap2idx1) %% 2 == 0)]
      hap2idx2 <- hap2idx2[which(1:length(hap2idx2) %% 2 == 0)]
      hap2 <- parenthap[,parent2idx[1]]
      for(j in 1:length(hap2idx1)){
        cropseg <- hap2idx1[j]:hap2idx2[j]
        cropseg <- cropseg[which(cropseg < m)]
        hap2[cropseg] <- parenthap[cropseg,parent2idx[2]]
      }
    }else{
      parent2idx <- which(colnames(parenthap) == indtbl$parent2[i])
      parent2idx <- sample(parent2idx, size=2, replace=F)
      hap2 <- parenthap[,parent2idx[1]]
    }
    return(c(hap1,hap2))
  }
  
  # Simulate individuals------------------------------------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  uniqchr <- unique(object$chr)
  chrsize <- rep(x = NA, times = length(uniqchr))
  names(chrsize) <- uniqchr
  nmkrchr <- chrsize
  for(i in 1:length(chrsize)){
    idx <- which(object$chr == uniqchr[i])
    bp <- object$bp[idx]
    chrsize[i] <- as.numeric(sum(diff(bp)))
    nmkrchr[i] <- length(idx)
  }
  idgen <- gsub(pattern = "( )|-|:", replacement = "", Sys.time())
  idgen <- paste0("ID",idgen,
                  sample(x = LETTERS, size = n.individuals, replace = TRUE),
                  sprintf(fmt = paste0("%0",as.integer(log10(n.individuals))+1,".f"), 1:n.individuals))
  indtbl <- data.frame(id = idgen,
                       parent1 = sample(x = names(parent1), size = n.individuals, prob = parent1, replace = TRUE),
                       parent2 = sample(x = names(parent2), size = n.individuals, prob = parent2, replace = TRUE),
                       stringsAsFactors = FALSE)
  write.table(x = indtbl, file = tmp.pedigree.file, col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
  write.table(x = cbind("SIM",idgen), file = tmp.samples.file, col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
  if("proportional" %in% model & length(model) == 1){
    chrprop <- chrsize/sum(chrsize)
    chrmean <- chrprop*length(chrsize)
    probs <- NULL
  }else if("uniform" %in% model & length(model) == 1){
    chrmean <- nmkrchr
    chrmean[1:length(chrmean)] <- 1
    probs <- NULL
  }else if(identical(names(model), names(chrsize)) & is.numeric(model) == TRUE){
    chrmean <- model*(chrsize/1e+6)/100
    model <- "chromosome"
    probs <- NULL
  }else if(sum(names(model) %in% object$marker) == length(model) & is.numeric(model) == TRUE){
    chrmean <- rep(x = NA, times = length(uniqchr))
    names(chrmean) <- uniqchr
    probs <- NULL
    for(k in 1:length(chrmean)){
      mkrtmp <- object$marker[which(object$chr == names(chrmean)[k])]
      mkrtmp <- mkrtmp[which(mkrtmp %in% names(model))]
      chrmean[k] <- mean(model[mkrtmp])*(chrsize[k]/1e+6)/100
      probs <- c(probs,model[mkrtmp]/sum(model[mkrtmp]))
    }
    model <- "marker"
  }else{
    emsg <- paste0("\nArgument model has to be one of the following:\n",
                   "a) 'uniform'\n",
                   "b) 'proportional'\n",
                   "c) named vector with chromosome-specific recombination rates\n",
                   "d) named vector with marker-specific recombination rates\n")
    stop(emsg)
  }
  for(chr in uniqchr){
    if(verbose == TRUE){
      cat("Simulating crossing over events on chromosome", names(nmkrchr[chr]), "\r")
    }
    m <- nmkrchr[chr]
    mkrsidx <- which(object$marker.in & object$chr == chr)
    mkrs <- object$marker[mkrsidx]
    if(model == "marker"){
      probchr <- probs[mkrs]
    }else{
      probchr <- NULL
    }
    parenthap <- ghap.slice(object = object, ids = unique(c(indtbl$parent1,indtbl$parent2)),
                            variants = mkrs, ncores = ncores)
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      inds <- parLapply(cl = cl, fun = indbuild, X = 1:n.individuals)
      stopCluster(cl)
    }else{
      inds <- mclapply(FUN = indbuild, X = 1:n.individuals, mc.cores = ncores)
    }
    inds <- matrix(data = unlist(inds), nrow = m, ncol = 2*n.individuals, byrow = F)
    mkrmap <- data.frame(CHR = object$chr[mkrsidx], MARKER = object$marker[mkrsidx], POS = object$bp[mkrsidx],
                         A0 = object$A0[mkrsidx], A1 = object$A1[mkrsidx], stringsAsFactors = FALSE)
    fwrite(x = as.data.table(mkrmap),
           file = tmp.markers.file, col.names = FALSE, row.names = FALSE,
           sep = " ", append = TRUE, nThread = ncores)
    fwrite(x = as.data.table(inds), file = tmp.phase.file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ",
           append = TRUE, nThread = ncores)
  }
  
  # Get files----------------------------------------------------------------------------------------
  if(verbose == TRUE){
    cat("Copying output files to the working directory... ")
  }
  ok <- file.copy(from = tmp.phase.file, to = phase.file)
  ok <- file.remove(tmp.phase.file)
  ok <- file.copy(from = tmp.markers.file, to = markers.file)
  ok <- file.remove(tmp.markers.file)
  ok <- file.copy(from = tmp.samples.file, to = samples.file)
  ok <- file.remove(tmp.samples.file)
  ok <- file.copy(from = tmp.pedigree.file, to = pedigree.file)
  ok <- file.remove(tmp.pedigree.file)
  if(verbose == TRUE){
    cat("Done.\n\n")
  }
  
}
