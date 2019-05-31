rm(list=ls(all=TRUE))
dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

library(mvtnorm)
library(reshape2)
library(MASS)
library(ggplot2)

### SQUARE EXPONENTIAL
K_SE <- function(x1,x2,l = 1) {exp(- abs(x1-x2)^2 / (l^2) ) }
### RATIONAL QUADRATIC (scale mixture of SE with different l)
K_RQ <- function(x1,x2,l= 1,alpha = 1) { ( 1 + abs(x1-x2)^2/(2 * alpha * l^2)  )^{-alpha} }
### OLS LINEAR REGRESSION
K_OLS <- function(x1,x2,sigma2 = 1){x1 * x2 + sigma2}
### PERIODIC WITH SE KERNEL
K_PR <- function(x1,x2,l = 1){exp(-2 * sin(pi*abs(x1 - x2))^2/(l^2))}


############ plotting
n = 20
t <- seq(0,10,len=n)
t_grid = expand.grid(t,t)
l_list = c(1,2,3,10,20,50)
nl = length(l_list)

function_Sigma_SE = function(l) matrix( apply(expand.grid(t,t),1, function(s) K_SE(x1=s[1],x2=s[2],l = l)),nrow = n)
Sigma_SE = lapply(l_list, function(l) function_Sigma_SE (l) )
Sigme_SE_melted = do.call('rbind',lapply(1:nl, function(i){df = melt(Sigma_SE[[i]]); df$parm = i;df$label ='SE';return(df) } ))
               

function_Sigma_RQ = function(l,a) matrix( apply(expand.grid(t,t),1, function(s) K_RQ(x1=s[1],x2=s[2],l = l, alpha = a)),nrow = n)
l_alpha_list = expand.grid(l_list[c(1,4,6)],c(0.1,2))
nla = nrow(l_alpha_list)
Sigma_RQ = lapply(1:nla, function(i) function_Sigma_RQ(l_alpha_list[i,1],l_alpha_list[i,2] ) )
Sigme_RQ_melted = do.call('rbind',lapply(1:nl, function(i){df = melt(Sigma_RQ[[i]]); df$parm = i;df$label ='RQ';return(df) } ))

function_Sigma_PR = function(l) matrix( apply(expand.grid(t,t),1, function(s) K_PR(x1=s[1],x2=s[2],l = l)),nrow = n)
Sigma_PR = lapply(l_list/5, function(l) function_Sigma_PR (l) )
Sigme_PR_melted = do.call('rbind',lapply(1:nl, function(i){df = melt(Sigma_PR[[i]]); df$parm = i;df$label ='PR';return(df) } ))


toplot_Sigma_st = do.call('rbind',list(Sigme_SE_melted,Sigme_RQ_melted,Sigme_PR_melted))
gg = ggplot(toplot_Sigma_st, aes(x = Var1, y = Var2,fill = value))
gg = gg + facet_grid(parm~label, scales = 'fixed' )
gg = gg + geom_tile()
gg = gg + scale_fill_gradientn(colors = c('lightyellow','blue','red'))
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'bottom')
gg = gg + labs(x = NULL, y = NULL)
ggsave(paste(dir,'/kernel_stationary.png',sep = ''),width =10,height = 8)

toplot_y_st_list = list()
for(s in 1:3){
  y_sim_RQ = do.call('rbind',lapply(1:nla, function(i){df = data.frame(t=t, y = mvrnorm(1, rep(0, n), Sigma_RQ[[i]] )); df$parm = i;df$label ='RQ';return(df) } ))
  y_sim_SE = do.call('rbind',lapply(1:nl, function(i){df = data.frame(t = t,y = mvrnorm(1, rep(0, n), Sigma_SE[[i]] )); df$parm = i;df$label ='SE';return(df) } ))
  y_sim_PR = do.call('rbind',lapply(1:nl, function(i){df = data.frame(t = t,y = mvrnorm(1, rep(0, n), Sigma_PR[[i]] )); df$parm = i;df$label ='PR';return(df) } ))
  toplot_y_st_list[[s]] =  do.call('rbind',list(y_sim_SE,y_sim_RQ,y_sim_PR))
  toplot_y_st_list[[s]]$s = s
}
toplot_y_st = do.call('rbind',toplot_y_st_list)

gg = ggplot(toplot_y_st, aes(x = t, y = y,color = factor(s)))
gg = gg + facet_grid(parm~label, scales = 'fixed' )
gg = gg + geom_line()
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'none')
gg = gg + labs(x = NULL, y = NULL)
ggsave(paste(dir,'/y_sim_stationary.png',sep = ''),width =10,height = 8)


id_obs = sample(1:n,4,replace= F)
t_obs = t[id_obs]
t_pred = t[-id_obs]
y_obs = subset(toplot_y_st, toplot_y_st$s == 1 & toplot_y_st$label=='RQ' & toplot_y_st$parm == 1 & toplot_y_st$t %in% t_obs)[,'y']

mu_pred_SE = 
mu_pred_RQ = 
mu_pred_PR = 
  
confidence_interval = sigma * 1.96
q_80 = qmvnorm(p = 0.95, tail = c("upper.tail"), sigma = function_Sigma_SE(1)[1:2,1:2])



#### WIENER PROCESS
K_WP <- function(x1,x2) min(x1,x2)


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

#### GIVEN THE 2D MAPPING TODO


