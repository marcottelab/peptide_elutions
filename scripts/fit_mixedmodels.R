library(tidyverse)
library(future)
library(furrr)

library(mixtools)
library(argparse)

parser <- ArgumentParser(description='Fit many mixed_models')
parser$add_argument('--elut_long', action='store',
    help='File of complete PSMs in tidy form (elut_long.csv)')
parser$add_argument('--outfile', dest='outfile', action='store', default = "mixed_models.csv",
    help='Outfile name')
parser$add_argument('--numgauss', action='store', type = 'integer', default=2,
    help='Number of gaussians to fit')

args = parser$parse_args()

mixmod <- function(df, numgauss, ...){
  d <- df$FractionID

  fit <- normalmixEM(d, k = numgauss) #try to fit two Gaussians

  data.frame(lambda = c(fit$lambda[1], fit$lambda[2]),
           mu = c(fit$mu[1], fit$mu[2]),
           sigma = c(fit$sigma[1], fit$sigma[2]),
           loglik = c(fit$loglik, fit$loglik),
           dist = c(1,2))


   out <- data.frame(#ProteinID = ProteinID,
            lambda1 = fit$lambda[1], 
            lambda2 = fit$lambda[2],
            mu1 = fit$mu[1], 
            mu2 = fit$mu[2],
            sigma1 = fit$sigma[1], 
            sigma2 = fit$sigma[2],
            loglik = fit$loglik)
   return(out)
}

safe_mixmod <- possibly(mixmod,otherwise = data.frame())


plan(multiprocess)
print(args$elut_long)


# only applied to elut_sel so far. 
mixed_models <-
  read_csv(args$elut_long) %>%
  filter(pepcount > 0) %>% #write_csv("elut_long.csv")
  mutate(peporder = Start) %>%
  filter(!is.na(peporder)) %>% #Explore why there are NA peporder in the first place
  split(list(.$ProteinID, .$ExperimentID)) %>%
    future_map_dfr(safe_mixmod, numgauss = args$numgauss, .id = "ProteinID.exp", .progress = TRUE) %>%
  separate(ProteinID.exp, into = c("ProteinID", "ExperimentID"), sep = "[.]") %>%
arrange(loglik)


mixed_models %>% write_csv(args$outfile)
