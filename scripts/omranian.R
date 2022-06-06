
#Segmentation of biological multivariate time series data
#Authors: Nooshin Omranian, Bernd Mueller-Roeber and and Zoran Nikoloski

#Sientific Reports


### Random data generation for segmentation problem



create_synthetic_data <- function()
{
  # how many clusters or groups per segment, i.e. 2 for first segment , 6 for 2nd etc...
  groupsPerInterval <- c(5,2,8,3,7,4)
  # length of the segment
  intervals <- c(5,8,12,3,7,5)
  # how many random variables for each group
  variablesEach <- 10
  
  test.data <- createRandomTimeSeries(intervals,groupsPerInterval,variablesEach,variableNoise=1)
  
  return(test.data)
}  

createRandomTimeSeries <- function(intervals,groupsPerInterval,variablesEach,variableNoise) {
  data <- matrix(0,nrow=max(groupsPerInterval),ncol=sum(intervals))
  for (i in 1:length(groupsPerInterval)) {
    end <- sum(intervals[1:i])
    if (i == 1) {
      start = 1
    } else {
      start = end - intervals[i] + 1
    }
    # print(start:end)
    tmp <- matrix(0,nrow=groupsPerInterval[i],ncol=length(start:end))
    for (j in 1:groupsPerInterval[i]) {
      tmp[j,] <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = length(start:end))[-1]
    }
    use <- getSample(1,groupsPerInterval[i],max(groupsPerInterval))
    # cat("...\n")
    # print(use)
    for (j in 1:length(use)) {
      data[j,start:end] <- tmp[use[j],]
    }
  }
  # smooth segment boundaries	
  data <- t(apply(data,1,function(x,leftBound) {
    for (i in leftBound) {
      x[i] <- (x[i-1]+2*x[i]+x[i+1])/4
    }
    x
  },cumsum(intervals)[-length(intervals)]))
  # add variables
  complete.data <- matrix(0,nrow=variablesEach*max(groupsPerInterval),ncol=sum(intervals))
  index <- 1
  for(i in 1:nrow(data)) {
    for (j in 1:variablesEach) {
      complete.data[index,] <- sapply(data[i,], function (x) {rnorm(1,mean=x,sd=variableNoise)})
      index <- index + 1
    }
  }
  return(complete.data)
}

getSample <- function(from,to,n) {
  notAll <- TRUE
  while(notAll) {
    this.sample <- sample(from:to,n,replace=TRUE)
    notAll <- !length(unique(this.sample)) == to
  }
  return(this.sample)
}



# Load the time series matrix including the comopnents profiles

# Here, the synthetic data is given as an example. To test the algorithm with synthetic data, uncomment the following two lines
# file <- "../Data/data_syn.obj"
# load(file)

MTS.Segmentation <- function(data, lambda1=NULL, lambda2=NULL)
{
  main.data <- data
  
  data <- data[,seq(from=dim(data)[2],to=1,by=-1)]
  
  C<-matrix(0,ncol(data),ncol(data))
  
  diag(C) <-0
  
  time <- t(data)
  
  i <- 1
  dd <- ncol(data)-2
  
  if (is.null(lambda1)==T | is.null(lambda1)==T)
  {
    lambda1 <- rep(0,dd)
    lambda2 <- rep(0,dd)
    
    while(i<= dd)
    {
      print(paste("step",i))
      
      y <- time[i,]
      
      x <- data[,-c(1)]
      
      a <- rep(1,ncol(x))
      
      if (ncol(x) <10)
        fld <- 5
      else
        fld <- 10
      if (ncol(x) <5)
        fld <- ncol(x)
      
      while(1)
      {
        opt1 <- NULL
        tryCatch(opt1 <- {withTimeout({optL1(response=y,penalized=x,fusedl=a,fold=fld,minlambda1=1,maxlambda1=50,trace=F,standardize=T);},
                                          timeout=3600,cpu=Inf)},
                 TimeoutException = function(ex) withTimeout({cat("Timeout. Skipping.\n")},timeout=Inf))
        if (is.null(opt1)==F)
          break
      }
      
      lambda1[i] <- opt1$lambda
      
      print(paste("lambda1",i, lambda1[i]))
      
      while(1)
      {
        opt2 <- NULL
        tryCatch(opt2 <- {withTimeout({optL2(response=y,penalized=x,fusedl=a,fold=fld,lambda1=opt1$lambda,minlambda2=1,maxlambda2=50,trace=F,standardize=T);},
                                          timeout=3600,cpu=Inf)},
                 TimeoutException = function(ex) withTimeout({cat("Timeout. Skipping.\n")},timeout=Inf))
        if (is.null(opt2)==F)
          break
      }
      
      lambda2[i] <- opt2$lambda
      
      print(paste("lambda2",i, lambda2[i]))
      
      fit <- penalized(y,x,lambda1=opt1$lambda,lambda2=opt2$lambda,fusedl=a,standardize=T,trace=F)
      
      t1 <- coefficients(fit,standardize=T,"all")[-1]
      
      C[i,(i+1):ncol(C)]<-t1
      
      data <- data[,-1]
      
      i<-i+1
    }
  }
  else
  {
    while(i<= dd)
    {
      y <- time[i,]
      
      x <- data[,-c(1)]
      
      a <- rep(1,ncol(x))
      
      if (ncol(x) <10)
        fld <- 5
      else
        fld <- 10
      if (ncol(x) <5)
        fld <- ncol(x)
      
      while(1)
      {
        fit <- NULL
        
        tryCatch(fit <- {withTimeout({penalized(y,x,lambda1=lambda1[i],lambda2=lambda2[i],fusedl=a,standardize=T,trace=F);},
                                         timeout=120,cpu=Inf)},
                 TimeoutException = function(ex) withTimeout({cat("Timeout. Skipping.\n")},timeout=Inf))  
        if (is.null(fit)==F)
          break 
      }
      
      t1 <- coefficients(fit,standardize=T,"all")[-1]
      
      C[i,(i+1):ncol(C)]<-t1
      
      data <- data[,-1]
      
      i<-i+1
    }
  }
  
  data <- main.data
  
  result <- list()
  
  C <- C[seq(from=dim(data)[2],to=1,by=-1),seq(from=dim(data)[2],to=1,by=-1)]
  
  result[[1]] <- C
  
  # neglecting the first and the last three time points
  if (ncol(data)<15)
    rm <- 2
  else
    rm <- 3
  C[1:rm,] <-0.01
  C[,1:rm] <-0
  C[,(dd+2-(rm-1)):(dd+2)] <-0
  C[(dd+2-(rm-1)):(dd+2),] <-0.01
  
  # calculating the absolute average sum of the triangular matrix C which inludes the regression coefficients
  A <- apply(abs(C),2,sum)
  A[c(1,dd+2)]=0
  r<-dd-(rm-1)*2
  A <- A/c(rep(r,rm),(r):1,rep(1,rm))
  
  local_minima <- which(diff(sign(diff(A)))==2)+1
  
  bps <- local_minima
  
  load(file)
  ts.plot(t(data),ylim=range(data),col="grey39")
  abline(v=bps,col="darkgreen",lwd=4)
  par(new=T)
  ts.plot(c(1),ylim=range(A),xlim=c(1,dd+2),gpars=list(xaxt="n",yaxt="n"),ylab="",xlab="")
  lines(A,type="o",col="red",lwd=2)
  
  result[[2]] <- bps
  
  result[[3]] <- lambda1
  
  result[[4]] <- lambda2
  
  return(result)
}