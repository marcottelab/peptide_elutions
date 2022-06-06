


wrap_locate_change<- function(df, noise_min = 0.01, noise_max = 0.5){
  
  df <- df %>% select(-peak)
  
  mat <- df %>% 
    data.matrix() %>% 
    addNoise(noise_min = noise_min, noise_max = noise_max) 
  
  change <- locate.change2(addNoise(mat, 0.01,0.5), standardize.series = TRUE, view.cusum = FALSE)
  change_tidy <- enframe(change) %>% pivot_wider(names_from = name, values_from = value) %>%
    unnest(c(changepoint, cusum)) %>%
    #mutate(peptide = names(a)[changepoint]) %>%
    mutate(leftr2 = leftfit(cusum.proj, changepoint),
           rightr2 = rightfit(cusum.proj, changepoint)) %>%
    mutate(pepidleft = names(df)[changepoint],
           pepidright = names(df)[changepoint + 1])
  
  return(change_tidy)
}

possibly_wrap_locate_change <- possibly(wrap_locate_change, otherwise = data.frame())


wrap_locate_change_repeat <- function(peak_sel){
  change_stats <- rerun(10,  possibly_wrap_locate_change(peak_sel))  %>%
    bind_rows(.id = 'run')
  return(change_stats)
}

# This is a modification of the one from  InspectBreakpoint that also returns the cumsum.proj
locate.change2 <- function (x, lambda, schatten = 2, sample.splitting = FALSE, 
                            standardize.series = FALSE, view.cusum = FALSE) 
{
  x <- as.matrix(x)
  if (dim(x)[2] == 1) 
    x <- t(x)
  p <- dim(x)[1]
  n <- dim(x)[2]
  if (missing(lambda)) 
    lambda <- sqrt(log(log(n) * p)/2)
  if (standardize.series) 
    x <- rescale.variance(x)
  if (sample.splitting) {
    x1 <- x[, seq(1, n, by = 2)]
    x2 <- x[, seq(2, n, by = 2)]
  }
  else {
    x1 <- x
    x2 <- x
  }
  cusum.matrix1 <- cusum.transform(x1)
  if (sample.splitting) {
    cusum.matrix2 <- cusum.transform(x2)
  }
  else {
    cusum.matrix2 <- cusum.matrix1
  }
  if (lambda >= max(abs(cusum.matrix1))) 
    lambda <- max(abs(cusum.matrix1)) - 1e-10
  vector.proj <- sparse.svd(cusum.matrix1, lambda, schatten)
  cusum.proj <- t(cusum.matrix2) %*% vector.proj
  if (view.cusum) 
    plot(as.numeric(cusum.proj), ylab = "projected cusum", 
         pch = 20)
  ret <- NULL
  ret$changepoint <- which.max(abs(cusum.proj))
  if (sample.splitting) 
    ret$changepoint <- ret$changepoint * 2
  ret$cusum <- max(abs(cusum.proj))
  ret$vector.proj <- vector.proj
  ret$cusum.proj <- as.numeric(cusum.proj)
  return(ret)
}


leftfit <- function(cusum.proj, changepoint){
  
  x <- rescale_0_1(cusum.proj[[1]])
  
  lefthalf <- x[1:changepoint]
  lm(seq(1:length(lefthalf)) ~ lefthalf) %>% summary() ->lh 
  leftmse <- rmse(lh)
  
  return(leftmse)
}

rmse <- function(sm) {
  sqrt(mean(sm$residuals^2))
}
rightfit <- function(cusum.proj, changepoint){
  
  x <- rescale_0_1(cusum.proj[[1]])
  righthalf <- x[changepoint:length(x)]
  lm(seq(1:length(righthalf)) ~ righthalf) %>% summary()->rh
  rightmse <- rmse(rh)
  return(rightmse)
}

rescale_0_1 <- function(x){
  (x-min(x))/(max(x)-min(x))
}


pwfit <- function(cusum.proj, changepoint){
  
  x <- rescale_0_1(cusum.proj[[1]])
  pw<- lm(seq(1:length(x)) ~ x*(x < changepoint) + x*(x>=changepoint)) %>% summary()
  #pwr <- pw$adj.r.squared 
  
  pwr <- rmse(pw)
  return(pwr)    
}


addNoise <- function(mat, noise_min = 0.01, noise_max = 1) {
  rand_mat <- matrix(runif(prod(dim(mat)), min = noise_min, max = noise_max), nrow = dim(mat)[1])
  return(mat + rand_mat)
  #https://stats.stackexchange.com/questions/46302/adding-noise-to-a-matrix-vector
  #if (!is.matrix(mtx)) mtx <- matrix(mtx, byrow = TRUE, nrow = 1)
  #random.stuff <- matrix(runif(prod(dim(mtx)), min = -(noise_amount), max = noise_amount), nrow = dim(mtx)[1])
  #random.stuff + mtx
}



# Get BayesBoundaries

get_boundaries <- function(samp){
  frac_vect <- samp %>% select(FractionOrder, pepcount) %>%
    uncount(pepcount) %>% pull(FractionOrder)
  gauss_intersections <-  AdaptGauss_nonshiny(frac_vect, fast = TRUE) 
  
  
  sse <- gauss_intersections$sse
  boundaries <- 
    tibble(FractionOrder = gauss_intersections$BayesBoundaries) %>% 
    rowid_to_column(var = "boundary") %>% 
    mutate(FractionOrder = as.integer(FractionOrder)) 
  
  onehundred_fractions <- tibble(FractionOrder = seq(100))
  #Can't join direction since breakpoint could be in a fraction without a peptide
  #print(onehundred_fractions)
  #print(boundaries)
  peaks <- onehundred_fractions %>%
    left_join(boundaries, by = "FractionOrder") %>% 
    mutate(peak = boundary) %>%
    fill(peak, .direction = "down") %>% 
    mutate(peak = replace_na(peak, 0 ))  %>%
    select(-boundary) %>%
    mutate(sse = sse)
  
  return(peaks)
}

possibly_get_boundaries <- possibly(get_boundaries, otherwise = data.frame())


getOptGauss_mod2 <- function(Data, Kernels, ParetoDensity,fast, maxgauss = 5){
  
  
  
  #Mean <- mean(Data)
  #Deviation <- sd(Data)
  #Weight <- 1
  Var_opt=EMGauss(Data,fast=fast)
  Means_opt <- Var_opt$Means
  Deviations_opt <- Var_opt$SDs
  Weights_opt <- Var_opt$Weights
  Fi_opt <- dnorm(Kernels,Means_opt,Deviations_opt)*Weights_opt
  RMS_opt <- sqrt(sum((Fi_opt-ParetoDensity)^2))
  
  ngauss <- seq(2,maxgauss,1)
  
  for(i in seq(1,10)){
      for(n in ngauss){
        Means <- rep(0 ,n)
        Deviations <-  rep(0 ,n)
        Weights <-  rep(0 ,n)
        Valskmeans <- kmeans(Data, n,iter.max=100)
        KValues <- Valskmeans$cluster
        
        for (i in 1:n){
          Means[i] <- mean(Data[KValues==i])
          Deviations[i] <- sd(Data[KValues==i])
          Weights[i] <- sum(KValues==i)
          if (is.na(Deviations[i])) {Deviations[i] <- 0}
        }
        Weights <- Weights/length(KValues)
        Var=EMGauss(Data,K=length(Means), Means,Deviations,Weights,10,fast=fast)
        Means <- Var$Means
        Deviations <- Var$SDs
        Weights <- Var$Weights
        Fi <- 0
        for (i in 1:n){
          Fi <- Fi+dnorm(Kernels,Means[i],Deviations[i])*Weights[i]
        }
        RMS <- sqrt(sum((Fi-ParetoDensity)^2))
        
        if(RMS < RMS_opt){
          print(n)
          
          Means_opt <- Means
          Deviations_opt <- Deviations
          Weights_opt <- Weights
          RMS_opt <- RMS
          print(RMS_opt)
          
        }
    
     }
    
  }
  
  order <- order(Means_opt)
  Means_opt <- Means_opt[order]
  Deviations_opt <- Deviations_opt[order]
  Weights_opt <- Weights_opt[order]
  out=list(means= Means_opt,deviations=Deviations_opt,weights=Weights_opt, rms = RMS_opt )
  return(out)
  
  

}




# Remove Shiny App parts of the AdaptGauss function from AdaptGauss package
# Modify Gaussian fitting to look for 1,2, or 3 Gaussians
AdaptGauss_nonshiny <- function (Data, Means = NaN, SDs = NaN, Weights = NaN, ParetoRadius = NaN, 
          LB = NaN, HB = NaN,  fast = T) 
{
  
                            data <- Data
                            GM <- Means
                            GS <- SDs
                            GW <- Weights
                            ParetoRadius <- ParetoRadius
                            LB <- LB
                            HB <- HB
                            nSignif <- 4
                            
                        
                            if (length(data) == 0) return("Error: Data could not be loaded")
                            
                            # CDM Remove Na/Nan  
                            dataNew <- 0
                            j <- 1
                            for (i in 1:length(data)) {
                              if (!is.nan(data[i]) && !is.na(data[i])) {
                                dataNew[j] <- data[i]
                                j <- j + 1
                              }
                            }
                            data <- dataNew
                            if (is.nan(LB)) LB <- min(data)
                            if (is.nan(HB)) HB <- max(data)
                            
                            # CDM Reduce to left and right bounds
                            dataNew <- 0
                            j <- 1
                            for (i in 1:length(data)) {
                              if (data[i] >= LB && data[i] <= HB) {
                                dataNew[j] <- data[i]
                                j <- j + 1
                              }
                            }
                            # CDM reduce to max sample size
                            data <- dataNew
                            if (length(data) > 25000) {
                              data <- data[rsampleAdaptGauss(25000, length(data))]
                              print("Reducing to 25000 datapoints")
                            }
                            if (is.nan(ParetoRadius)) {
                              ParetoRadius <- DataVisualizations::ParetoRadius(data)
                              nRow = length(data)
                            }
                            ParetoDensityEstimationVar <- DataVisualizations::ParetoDensityEstimation(data, 
                                                                                                      paretoRadius = ParetoRadius)
                            ParetoDensity <- ParetoDensityEstimationVar$paretoDensity
                            Kernels <- ParetoDensityEstimationVar$kernels
                           
                            
                            # CDM get the gaussians (1,2,3, or 4)
                            if (is.nan(sum(GM)) || is.nan(sum(GS)) || is.nan(sum(GW)) || 
                                length(GM) != length(GS) || length(GM) != length(GW)) {
                              vars = getOptGauss_mod2(Data = data, Kernels, ParetoDensity, 
                                                 fast = fast, maxgauss = 5)
                              GM <- vars$means
                              GS <- vars$deviations
                              GW <- vars$weights
                              rms <- vars$rms
                            }
                            
                            BB <- NaN
                            numGauss <- length(GM)
                            
                            #meanRMS0 <- mean(data)
                            #DeviationRMS0 <- sd(data)
                            #Fi <- dnorm(Kernels, meanRMS0, DeviationRMS0)
                            #RMS0 <- sqrt(sum((Fi - ParetoDensity)^2))
                            
                         
                            #RMS <- RMS0 # Added CDM
                              if (numGauss > 1 && sum(GW > 0) > 1) {
                                BayesBoundaries <- BayesDecisionBoundaries(GM[GW > 
                                                                                0], GS[GW > 0], GW[GW > 0])
                                BB <- BayesBoundaries # Changed from superassignment
                               
                              output <- list(Means = GM, SDs = GS, Weights = GW, 
                                             ParetoRadius = ParetoRadius, BayesBoundaries = BB, rms = rms)
                            
                         
                              } else if(numGauss ==1){
                                output <-list(Means = GM, SDs = GS, Weights = GW, rms = rms)
                              }
                           
                         
  return(output)
}


# Modified to test for 1,2, 3, or 4 gaussians
# Way too complicated, modify
getOptGauss_mod <- function(Data, Kernels, ParetoDensity,fast){ 
  # Teste RMS fuer einen Gauss
  Mean1 <- mean(Data)
  Deviation1 <- sd(Data)
  Weight1 <- 1
  Var=EMGauss(Data,fast=fast)
  Mean1 <- Var$Means
  Deviation1 <- Var$SDs
  Weight1 <- Var$Weights
  Fi <- dnorm(Kernels,Mean1,Deviation1)*Weight1
  RMS1 <- sqrt(sum((Fi-ParetoDensity)^2))
  
  # CDM Teste RMS fuer 2 Gauss
  
  Means2 <- c(0,0)
  Deviations2 <- c(0,0)
  Weights2 <- c(0,0)
  Valskmeans <- kmeans(Data,2,iter.max=1000)
  KValues <- Valskmeans$cluster
  
  
  # So pareto density comes in (fixed, and calculate volume of each gaussian and subtract)
  for (i in 1:2){
    Means2[i] <- mean(Data[KValues==i])
    Deviations2[i] <- sd(Data[KValues==i])
    Weights2[i] <- sum(KValues==i)
    if (is.na(Deviations2[i])) {Deviations2[i] <- 0}
  }
  Weights2 <- Weights2/length(KValues)
  Var=EMGauss(Data,K=length(Means2), Means2,Deviations2,Weights2,10,fast=fast)
  Means2 <- Var$Means
  Deviations2 <- Var$SDs
  Weights2 <- Var$Weights
  Fi <- 0
  for (i in 1:2){
    Fi <- Fi+dnorm(Kernels,Means2[i],Deviations2[i])*Weights2[i]
  }
  RMS2 <- sqrt(sum((Fi-ParetoDensity)^2))
  
  
  # Teste RMS fuer 3 Gauss
  Means3 <- c(0,0,0)
  Deviations3 <- c(0,0,0)
  Weights3 <- c(0,0,0)
  Valskmeans <- kmeans(Data,3,iter.max=100)
  KValues <- Valskmeans$cluster
  
  for (i in 1:3){
    Means3[i] <- mean(Data[KValues==i])
    Deviations3[i] <- sd(Data[KValues==i])
    Weights3[i] <- sum(KValues==i)
    if (is.na(Deviations3[i])) {Deviations3[i] <- 0}
  }
  Weights3 <- Weights3/length(KValues)
  Var=EMGauss(Data,K=length(Means3), Means3,Deviations3,Weights3,10,fast=fast)
  Means3 <- Var$Means
  Deviations3 <- Var$SDs
  Weights3 <- Var$Weights
  Fi <- 0
  for (i in 1:3){
    Fi <- Fi+dnorm(Kernels,Means3[i],Deviations3[i])*Weights3[i]
  }
  RMS3 <- sqrt(sum((Fi-ParetoDensity)^2))
  
  # Test for 4 gauss
  # Teste RMS fuer 3 Gauss
  Means4 <- c(0,0,0,0)
  Deviations4 <- c(0,0,0,0)
  Weights4 <- c(0,0,0,0)
  Valskmeans <- kmeans(Data,4,iter.max=100)
  KValues <- Valskmeans$cluster
  
  for (i in 1:4){
    Means4[i] <- mean(Data[KValues==i])
    Deviations4[i] <- sd(Data[KValues==i])
    Weights4[i] <- sum(KValues==i)
    if (is.na(Deviations4[i])) {Deviations4[i] <- 0}
  }
  Weights4 <- Weights4/length(KValues)
  Var=EMGauss(Data,K=length(Means4), Means4,Deviations4,Weights4,10,fast=fast)
  Means4 <- Var$Means
  Deviations4 <- Var$SDs
  Weights4 <- Var$Weights
  Fi <- 0
  for (i in 1:4){
    Fi <- Fi+dnorm(Kernels,Means4[i],Deviations4[i])*Weights4[i]
  }
  RMS4 <- sqrt(sum((Fi-ParetoDensity)^2))
  
  
  
  # Test for 4 gauss
  # Teste RMS fuer 3 Gauss
  Means5 <- c(0,0,0,0 ,0)
  Deviations5 <- c(0,0,0,0,0)
  Weights5 <- c(0,0,0,0,0)
  Valskmeans <- kmeans(Data,5,iter.max=100)
  KValues <- Valskmeans$cluster
  
  for (i in 1:5){
    Means5[i] <- mean(Data[KValues==i])
    Deviations5[i] <- sd(Data[KValues==i])
    Weights5[i] <- sum(KValues==i)
    if (is.na(Deviations5[i])) {Deviations5[i] <- 0}
  }
  Weights5 <- Weights5/length(KValues)
  Var=EMGauss(Data,K=length(Means5), Means5,Deviations5,Weights5,10,fast=fast)
  Means5 <- Var$Means
  Deviations5 <- Var$SDs
  Weights5 <- Var$Weights
  Fi <- 0
  for (i in 1:5){
    Fi <- Fi+dnorm(Kernels,Means5[i],Deviations5[i])*Weights5[i]
  }
  RMS5 <- sqrt(sum((Fi-ParetoDensity)^2))
  
  
  
  # ueberpruefe ob RMS1( 1 Gauss) oder RMS2 (3 Gauss ) kleiner ist. Speichere zugehoerige means, deviations und weights
  # check whether RMS1 (1 Gauss), RMS2 (2 Gauss) or RMS3 (3 Gauss) is smaller. Save associated means, deviations and weights
  # 3 is the number of parameters( mean, sd, weight)
  # Mean square errors * log(#parameters * #models)
  
  
  SSE1=RMS1^2*log(3)
  SSE2=RMS2^2*log(3*2)
  SSE3=RMS3^2*log(3*3)
  SSE4=RMS4^2*log(3*4)
  SSE5=RMS5^2*log(3*5)
  print(paste(RMS1, SSE1))
  print(paste(RMS2,SSE2))
  print(paste(RMS3, SSE3))
  print(paste(RMS4, SSE4))
  print(paste(RMS5, SSE5))
  if (RMS1 == min(RMS1, RMS2, RMS3, RMS4, RMS5)){ 
    #if(SSE1 == min(SSE1, SSE2, SSE3, SSE4)){
    print("1wins")
    means <- Mean1
    deviations <- Deviation1
    weights <- Weight1
    rms <- RMS1
    
  } else if(RMS2 == min(RMS1, RMS2, RMS3, RMS4, RMS5)){
    #} else if(SSE2 == min(SSE1, SSE2, SSE3, SSE4)){
    print("2wins")
    means <- Means2
    deviations <- Deviations2
    weights <- Weights2
    rms <- RMS2
  } else if(RMS3 == min(RMS1, RMS2, RMS3, RMS4, RMS5)){
    #} else if(SSE3 == min(SSE1, SSE2, SSE3, SSE4)){
    print("3wins")
    means <- Means3
    deviations <- Deviations3
    weights <- Weights3
    rms <- RMS3
  } else if(RMS4 == min(RMS1, RMS2, RMS3, RMS4, RMS5)){
    #} else if(SSE4 == min(SSE1, SSE2, SSE3, SSE4)){
    print("4wins")
    means <- Means4
    deviations <- Deviations4
    weights <- Weights4
    rms <- RMS4
  }
  else if(RMS5 == min(RMS1, RMS2, RMS3, RMS4, RMS5)){
    #} else if(SSE4 == min(SSE1, SSE2, SSE3, SSE4)){
    print("5wins")
    means <- Means5
    deviations <- Deviations5
    weights <- Weights5
    rms <- RMS5
  }
  # Ordne gaussians nach mean
  order <- order(means)
  means <- means[order]
  deviations <- deviations[order]
  weights <- weights[order]
  out=list(means=means,deviations=deviations,weights=weights, rms = rms )
  
  return(out)
}

