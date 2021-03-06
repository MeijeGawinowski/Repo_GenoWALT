rm(list=ls())

NB_Ne <- function (infile, alleles, sample.interval, bound = c(50, 1e+07), 
                   profile.likelihood = FALSE) {
  if (length(sample.interval) < 2) {
    stop("sample.interval must be a vector of length 2 or more. See help for details.")
  }
  if (any(sample.interval < 0)) {
    stop("sample.interval must be a non-negative vector. See help for details.")
  }
  if (any(sample.interval%%1 != 0)) {
    stop("sample.interval must contain non-negative integers. See help for details.")
  }
  if (any(alleles%%1 != 0)) {
    stop("alleles must be positive integers. See help for details.")
  }
  dirmulti.infile <- function(infile, alleles, sample.interval) {
    f <- function(dat.line) {
      temp <- strsplit(dat.line, " ")[[1]]
      temp <- as.numeric(temp[temp != ""])
      return(temp)
    }
    dat <- readLines(infile)
    dat <- lapply(dat, f)
    dat.length <- sapply(dat, length)
    dat <- dat[dat.length != 0]
    id <- (1:length(dat))%%length(alleles)
    temporal.samples <- length(sample.interval)
    if (length(dat) != temporal.samples * length(alleles)) {
      stop("infile format. See help for details.")
    }
    output <- list()
    for (i in 1:length(alleles)) {
      output[[i]] <- matrix(unlist(dat[id == (i - 1)]), 
                            ncol = alleles[i], byrow = T)
    }
    return(output)
  }
  ddirmulti <- function(x, alpha, log = T) {
    if (any(x%%1 != 0)) {
      stop("allele counts must be non-negative integers")
    }
    if (any(x < 0)) {
      stop("allele counts must be non-negative integers")
    }
    temp <- lgamma(sum(x) + 1) - sum(lgamma(x + 1)) + lgamma(sum(alpha)) - 
      lgamma(sum(alpha) + sum(x)) + sum(lgamma(alpha + 
                                                 x) - lgamma(alpha))
    if (log == T) {
      return(temp)
    }
    else {
      return(exp(temp))
    }
  }
  dirichlet.updating <- function(dat.list, N.dip, sample.interval = sample.interval) {
    dirichlet.parms <- dat.list * 0
    dirichlet.parms[1, ] <- 1
    for (i in 2:temporal.samples) {
      time.diff <- sample.interval[i] - sample.interval[i - 
                                                          1]
      kt <- (1 - 1/(2 * N.dip))^time.diff
      drift.parms <- kt/(1 - kt)
      temp <- dirichlet.parms[i - 1, ] + dat.list[i - 1, 
                                                  ]
      dirichlet.parms[i, ] <- temp * drift.parms/(1 + sum(temp) + 
                                                    drift.parms)
    }
    return(dirichlet.parms)
  }
  dirichlet.log.likelihood <- function(dat.list, N.dip, sample.interval) {
    dirichlet.parms <- dirichlet.updating(dat.list = dat.list, 
                                          N.dip = N.dip, sample.interval = sample.interval)
    likelihood.value <- rep(0, nrow(dat.list))
    for (i in 2:temporal.samples) {
      likelihood.value[i] <- ddirmulti(x = dat.list[i, 
                                                    ], alpha = dirichlet.parms[i, ])
    }
    return(sum(likelihood.value))
  }
  lapply.wrapper <- function(N.dip, dat, z = 0) {
    log.likelihood.overall <- lapply(dat, dirichlet.log.likelihood, 
                                     N.dip = N.dip, sample.interval)
    return(sum(unlist(log.likelihood.overall)) - z)
  }
  temporal.samples <- length(sample.interval)
  dat <- dirmulti.infile(infile = infile, alleles = alleles, 
                         sample.interval = sample.interval)
  result <- optimize(lapply.wrapper, interval = bound, dat = dat, 
                     maximum = T, tol = .Machine$double.eps^0.1)
  N.point <- result$maximum
  log.likelihood <- result$objective
  N.lb <- min(bound)
  try(N.lb <- uniroot(lapply.wrapper, c(min(bound), N.point), 
                      dat = dat, z = log.likelihood - 2, tol = .Machine$double.eps^0.15)$root, 
      silent = T)
  N.ub <- max(bound)
  try(N.ub <- uniroot(lapply.wrapper, c(N.point, max(bound)), 
                      dat = dat, z = log.likelihood - 2, tol = .Machine$double.eps^0.15)$root, 
      silent = T)
  if (profile.likelihood == TRUE) {
    N.value <- seq(N.lb, N.ub, length.out = 100)
    profile.CI <- cbind(N.value, sapply(N.value, lapply.wrapper, 
                                        dat = dat))
    colnames(profile.CI) <- c("log.like", "N")
    return(list(N = N.point, CI = c(N.lb, N.ub), log.like = log.likelihood, 
                profile.CI = profile.CI))
  }
  return(list(N = N.point, CI = c(N.lb, N.ub), log.like = log.likelihood))
}


N_rep <- 1000
v_Ne <- c()
for (rep in seq(0,N_rep-1)){
  print(rep)
  path_rep <- paste("~/Documents/GenoWALT/V4_tests/Tests/NE/LOC5/Fitness/Results/sizeA/Rep",rep,sep="")
  setwd(path_rep)
  res=data.frame() # tab for Rep 0
  N_loc=5
  N_gen=10
  for (gen in seq(0,N_gen)){
    file=paste("pop",gen,".txt",sep="")
    dat=read.table(file,header=TRUE)
    for (i in seq(1,N_loc*2,2)){
      v <- c(dat[,i],dat[,i+1]) # alleles at locus i
      n0i0 <- length(which(v==0)) # frequence allele 0 locus i time 0
      n1i0 <- length(which(v==1)) # frequence allele 1 locus i time 0
      row <- c(n0i0,n1i0)
      res <- rbind(res,row)
    }
  }
  
  write.table(res,"res.txt",sep=" ",col.names=FALSE,row.names=FALSE)
  
  alleles_res=rep(2,N_loc)
  int_res=seq(0,N_gen)
  
  res_NB <- NB_Ne(infile="res.txt",alleles=alleles_res,sample.interval=int_res)
  v_Ne <- c(v_Ne,res_NB$N)
}
