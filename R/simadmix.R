#Function: ghap.simadmix
#License: GPLv3 or later
#Modification date: 3 Jun 2022
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Simulate admixed individuals

ghap.simadmix <- function(
  object,
  n.individuals,
  n.generations,
  ancestors,
  proportions = NULL,
  alpha = NULL,
  out.file,
  only.active.markers = TRUE,
  ncores = 1,
  verbose = TRUE
){
  
  # Check if phase is a GHap.phase object --------------------------------------
  if(inherits(object, "GHap.phase") == FALSE){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  # Check if inactive markers should be reactivated ----------------------------
  if(only.active.markers == FALSE){
    object$marker.in <- rep(TRUE,times=object$nmarkers)
    object$nmarkers.in <- length(which(object$marker.in))
  }
  
  # Check if ancestors exist ---------------------------------------------------
  ancestors <- ancestors[sort(names(ancestors))]
  n.ancestors <- length(unlist(ancestors))
  K <- length(ancestors)
  if(length(which(unlist(ancestors) %in% object$id)) != n.ancestors){
    stop("\nSome of the ancestors are not present in the GHap.phase object.")
  }
  
  # Check if out file exist ----------------------------------------------------
  samples.file <- paste(out.file,"samples",sep=".")
  phase.file <- paste(out.file,"phase",sep=".")
  markers.file <- paste(out.file,"markers",sep=".")
  prop.file <- paste(out.file,"proportions",sep=".")
  haplo.file <- paste(out.file,"haplotypes",sep=".")
  if(file.exists(samples.file) == TRUE | file.exists(markers.file) == TRUE | file.exists(phase.file) == TRUE){
    stop("Output file already exists!")
  }else{
    rnumb <- runif(n = 1, min = 1, max = 1e+6)
    rnumb <- ceiling(rnumb)
    tmp.file <- paste(tempdir(),"/tmp",rnumb,sep="")
    tmp.samples.file <- paste(tmp.file,"samples",sep=".")
    tmp.phase.file <- paste(tmp.file,"phase",sep=".")
    tmp.markers.file <- paste(tmp.file,"markers",sep=".")
    tmp.prop.file <- paste(tmp.file,"proportions",sep=".")
    tmp.haplo.file <- paste(tmp.file,"haplotypes",sep=".")
  }
  
  # Check for ancestry proportions ---------------------------------------------
  rdir <- function(n, a){
    npar <- length(a)
    x <- matrix(data = rgamma(n = npar*n, shape = a), ncol = npar, byrow = TRUE)
    s <- as.numeric(x%*%rep(x = 1, times = npar))
    x <- x/s
    return(x)
  }
  if(is.null(alpha) == FALSE & is.null(proportions) == FALSE){
    stop("\nProvide alpha values OR a proportions data frame - not both.")
  }
  if(is.null(alpha) == FALSE & is.null(proportions) == TRUE){
    if(identical(sort(names(alpha)), names(ancestors)) == FALSE){
      stop("\nThe list of alpha values do not match ancestors.")
    }
    alpha <- alpha[names(ancestors)]
    proportions <- rdir(n = n.individuals, a = unlist(alpha))
    proportions <- as.data.frame(proportions)
    colnames(proportions) <- names(ancestors)
    method <- "sampled from a Dirichelet distribution"
  }
  if(is.null(alpha) == TRUE & is.null(proportions) == FALSE){
    if(nrow(proportions) != n.individuals){
      stop("\nNumber of rows in proportions do not match the number of individuals.")
    }
    if(ncol(proportions) != K){
      stop("\nNumber of columns in proportions do not match the number of ancestral populations.")
    }
    proportions <- proportions[,names(ancestors)]
    method <- "provided by the user"
  }
  if(is.null(alpha) == TRUE & is.null(proportions) == TRUE){
    proportions <- rdir(n = n.individuals, a = rep(x = 1, times = K))
    proportions <- as.data.frame(proportions)
    colnames(proportions) <- names(ancestors)
    method <- "sampled from a Dirichelet distribution"
  }
  
  # Settings of ancestry simulation --------------------------------------------
  if(verbose == TRUE){
    cat("\n\nAncestry simulation started...\n")
    cat("Number of individuals: ", n.individuals, "\n", sep="")
    cat("Number of ancestors: ", n.ancestors, "\n", sep="")
    cat("Number of ancestral populations: ", K, "\n", sep="")
    cat("Number of generations since admixture: ", n.generations, "\n", sep="")
    cat("Ancestry proportions ", method, "\n\n", sep="")
  }
  
  # Function for individual simulation -----------------------------------------
  indbuild <- function(id){
    newhap <- list(NULL,NULL)
    ancsmooth <- NULL
    for(k in 1:2){
      j <- 0
      while(j < m){
        i = j + 1
        p <- sample(1:K, size = 1, prob = proportions[id,])
        h <- sample(x = unlist(ancestors[p]), size = 1)
        h <- which(colnames(anchap) == h)
        h <- sample(x = h, size = 1)
        parent <- anchap[,h]
        x <- ceiling(rexp(n = 1, rate = 1/l))
        j <- min(c(i + x, m))
        frag <- paste(parent[i:j], collapse = "")
        newhap[[k]] <- paste(newhap[[k]],frag,sep="")
        seg <- paste(c("SIM",indtbl$ID[id],k,chr,bp[i],bp[j],bp[j]-bp[i]+1,names(ancestors)[p]),collapse="\t")
        ancsmooth <- paste(ancsmooth,seg,sep="\n")
      }
    }
    ancsmooth <- gsub(pattern = "^\n", replacement = "", x = ancsmooth)
    return(c(ancsmooth,newhap[[1]],newhap[[2]]))
  }
  indformat <- function(id){
    hap1 <- as.numeric(unlist(strsplit(x = inds[id,2], split = "")))
    hap2 <- as.numeric(unlist(strsplit(x = inds[id,3], split = "")))
    return(c(hap1,hap2))
  }
  
  # Simulate individuals -------------------------------------------------------
  ncores <- min(c(detectCores(), ncores))
  uniqchr <- unique(object$chr)
  chrsize <- rep(x = NA, times = length(uniqchr))
  names(chrsize) <- uniqchr
  nmkrchr <- chrsize
  chrdens <- chrsize
  for(i in 1:length(chrsize)){
    idx <- which(object$chr == uniqchr[i] & object$marker.in)
    bp <- object$bp[idx]
    mrkdist <- diff(bp)
    chrsize[i] <- as.numeric(sum(mrkdist))
    nmkrchr[i] <- length(idx)
    chrdens[i] <- 1e+6/mean(mrkdist)
  }
  idgen <- gsub(pattern = "( )|-|:", replacement = "", Sys.time())
  idgen <- paste0("ID",idgen,
                  sample(x = LETTERS, size = n.individuals, replace = TRUE),
                  sprintf(fmt = paste0("%0",as.integer(log10(n.individuals))+1,".f"), 1:n.individuals))
  idgen <- sort(idgen)
  indtbl <- cbind("SIM", idgen, proportions)
  colnames(indtbl)[1:2] <- c("POP","ID")
  write.table(x = cbind("SIM",idgen), file = tmp.samples.file, col.names = FALSE,
              row.names = FALSE, sep = " ", quote = FALSE)
  write.table(x = indtbl, file = tmp.prop.file, col.names = TRUE,
              row.names = FALSE, sep = " ", quote = FALSE)
  write(x = "POP\tID\tHAP\tCHR\tBP1\tBP2\tSIZE\tANCESTRY", file = tmp.haplo.file)
  for(chr in uniqchr){
    if(verbose == TRUE){
      cat("Simulating admixture on chromosome", names(nmkrchr[chr]), "\r")
    }
    m <- nmkrchr[chr]
    mkrsidx <- which(object$marker.in & object$chr == chr)
    mkrs <- object$marker[mkrsidx]
    anchap <- ghap.slice(object = object, ids = unlist(ancestors),
                         variants = mkrs, ncores = ncores)
    l <- (100/(2*n.generations))*chrdens[chr]
    if(Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(ncores)
      inds <- parLapply(cl = cl, fun = indbuild, X = 1:n.individuals)
      stopCluster(cl)
      inds <- matrix(data = unlist(inds), ncol = 3, byrow = TRUE)
      anctracks <- inds[,1]
      cl <- makeCluster(ncores)
      inds <- parLapply(cl = cl, fun = indformat, X = 1:n.individuals)
      stopCluster(cl)
      inds <- matrix(data = unlist(inds), nrow = m, ncol = 2*n.individuals, byrow = FALSE)
    }else{
      inds <- mclapply(FUN = indbuild, X = 1:n.individuals, mc.cores = ncores)
      inds <- matrix(data = unlist(inds), ncol = 3, byrow = TRUE)
      anctracks <- inds[,1]
      inds <- mclapply(FUN = indformat, X = 1:n.individuals, mc.cores = ncores)
      inds <- matrix(data = unlist(inds), nrow = m, ncol = 2*n.individuals, byrow = FALSE)
    }
    mkrmap <- data.frame(CHR = object$chr[mkrsidx], MARKER = object$marker[mkrsidx],
                         POS = object$bp[mkrsidx], A0 = object$A0[mkrsidx],
                         A1 = object$A1[mkrsidx], stringsAsFactors = FALSE)
    fwrite(x = as.data.table(mkrmap),
           file = tmp.markers.file, col.names = FALSE, row.names = FALSE,
           sep = " ", append = TRUE, nThread = ncores)
    fwrite(x = as.data.table(inds), file = tmp.phase.file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ",
           append = TRUE, nThread = ncores)
    write(x = anctracks, file = tmp.haplo.file, append = TRUE)
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
  ok <- file.copy(from = tmp.prop.file, to = prop.file)
  ok <- file.remove(tmp.prop.file)
  ok <- file.copy(from = tmp.haplo.file, to = haplo.file)
  ok <- file.remove(tmp.haplo.file)
  if(verbose == TRUE){
    cat("Done.\n\n")
  }
  
}
