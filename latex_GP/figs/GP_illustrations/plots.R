rm(list=ls(all=TRUE))
dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

library(mvtnorm)
library(reshape2)

### SQUARE EXPONENTIAL
K_SE <- function(x1,x2,l = 1) {exp(- abs(x1-x2)^2 / (l^2) ) }
### RATIONAL QUADRATIC (scale mixture of SE with different l)
K_RQ <- function(x1,x2,l= 1,alpha = 1) { ( 1 + abs(x1-x2)^2/(2 * alpha * l^2)  )^{-alpha} }
### OLS LINEAR REGRESSION
K_OLS <- function(x1,x2,sigma2 = 1){x1 * x2 + sigma2}
### PERIODIC WITH SE KERNEL
K_PR <- function(x1,x2,l = 1){exp(-2 * sin(pi*abs(x1 - x2))^2/(l^2))}

#### NN with one hidden layer (non stationary)

K_NN <- function(x1,x2,Sigma){ 
  A = matrix(c(1,x1),ncol = 1) %*% Sigma %*% matrix(c(1,x2),nrow = 1)
  B = 1 + matrix(c(1,x1),ncol = 1) %*% Sigma %*% matrix(c(1,x1),nrow = 1)
  C = 1 + matrix(c(1,x2),ncol = 1) %*% Sigma %*% matrix(c(1,x2),nrow = 1)
  return( 2 * asin(2 * A / sqrt(B * C) ) / pi )
}
#### GIVEN GAUSSIAN BASIS FUNCTION
K_G <- function(x1,x2,sigma2_u, sigma2_g){
  sigma2_e_inv = 2/sigma2_g + 1/sigma2_u
  sigma2_s = 2 * sigma2_g + sigma2_g^2/sigma2_u
  sigma2_m = 2 * sigma2_u + sigma2_g
  A = exp( - x1*x1/(2* sigma2_m))
  B = exp( - x2*x2/(2* sigma2_m))
  C = exp( - (x1 - x2)^2/(2* sigma2_s))
  return(1/(sigma2_e_inv * sigma2_u)  * A * B * C)
}

#### GIVEN THE 2D MAPPING TO  

############ plotting
n = 20
t <- seq(0,10,len=n)
t_grid = expand.grid(t,t)
l_list = c(1,2,3,10,20,100)
nl = length(l_list)

function_Sigma_SE = function(l) matrix( apply(expand.grid(t,t),1, function(s) K_SE(x1=s[1],x2=s[2],l = l)),nrow = n)

Sigma_SE = lapply(l_list, function(l) function_Sigma_SE (l) )
Sigme_SE_melted = do.call('rbind',lapply(1:nl, function(i){df = melt(Sigma_SE[[i]]); df$parm = 1;df$label ='SE';return(df) } ))




function_Sigma_SE(1)


q_80 = qmvnorm(p = 0.8, tail = c("upper.tail"), sigma = function_Sigma_SE(1))

