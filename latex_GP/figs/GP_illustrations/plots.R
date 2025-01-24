rm(list=ls(all=TRUE))
dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

library(mvtnorm)
library(reshape2)
library(MASS)
library(ggplot2)

### SQUARE EXPONENTIAL
K_SE <- function(x1,x2,l = 1) {exp(- abs(x1-x2)^2 / (2 * l^2) ) }
### RATIONAL QUADRATIC (scale mixture of SE with different l)
K_RQ <- function(x1,x2,l= 1,alpha = 1) { ( 1 + abs(x1-x2)^2/(2 * alpha * l^2)  )^{-alpha} }
### OLS LINEAR REGRESSION
K_OLS <- function(x1,x2,sigma2 = 1){x1 * x2 + sigma2}
### PERIODIC WITH SE KERNEL
K_PR <- function(x1,x2,l = 1,p = 1){exp(-2 * sin(pi*abs(x1 - x2) / p)^2 / (l^2) ) }
### LOCALLY PERIODIC
K_locPR <- function(x1,x2,l = 1,p = 1) {K_PR(x1,x2,l,p) * K_SE(x1,x2,l)}


############ plotting
n = 20
t <- seq(0,10,len=n)
t_grid = expand.grid(t,t)

l_list = c(0.25,1,3,0.25,1,3)
nl = length(l_list)

l_alpha_list = expand.grid(l_list[c(1,2,3)],c(0.5,7))
nla = nrow(l_alpha_list)

function_Sigma_SE = function(l) matrix( apply(expand.grid(t,t),1, function(s) K_SE(x1=s[1],x2=s[2],l = l)),nrow = n)
Sigma_SE = lapply(l_list, function(l) function_Sigma_SE (l) )
Sigme_SE_melted = do.call('rbind',lapply(1:nl, function(i){df = melt(Sigma_SE[[i]]); df$parm = i;df$label ='SE';return(df) } ))
               
function_Sigma_RQ = function(l,a) matrix( apply(expand.grid(t,t),1, function(s) K_RQ(x1=s[1],x2=s[2],l = l, alpha = a)),nrow = n)
Sigma_RQ = lapply(1:nla, function(i) function_Sigma_RQ(l_alpha_list[i,1],l_alpha_list[i,2] ) )
Sigme_RQ_melted = do.call('rbind',lapply(1:nl, function(i){df = melt(Sigma_RQ[[i]]); df$parm = i;df$label ='RQ';return(df) } ))

function_Sigma_PR = function(l,a) matrix( apply(expand.grid(t,t),1, function(s) K_PR(x1=s[1],x2=s[2],l = l, p = a)),nrow = n)
Sigma_PR = lapply(1:nla, function(i) function_Sigma_PR(l_alpha_list[i,1],l_alpha_list[i,2] /2) )
Sigme_PR_melted = do.call('rbind',lapply(1:nl, function(i){df = melt(Sigma_PR[[i]]); df$parm = i;df$label ='PR';return(df) } ))

function_Sigma_locPR = function(l,a) matrix( apply(expand.grid(t,t),1, function(s) K_locPR(x1=s[1],x2=s[2],l = l, p = a)),nrow = n)
Sigma_locPR = lapply(1:nla, function(i) function_Sigma_locPR(l_alpha_list[i,1],l_alpha_list[i,2]/2 ) )
Sigme_locPR_melted = do.call('rbind',lapply(1:nl, function(i){df = melt(Sigma_locPR[[i]]); df$parm = i;df$label ='locPR';return(df) } ))


toplot_Sigma_st = do.call('rbind',list(Sigme_SE_melted,Sigme_RQ_melted,Sigme_PR_melted,Sigme_locPR_melted))
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
ggsave(paste(dir,'/kernel_stationary.eps',sep = ''),width =10,height = 8)

toplot_y_st_list = list()
for(s in 1:3){
  y_sim_RQ = do.call('rbind',lapply(1:nla, function(i){df = data.frame(t=t, y = mvrnorm(1, rep(0, n), Sigma_RQ[[i]] )); df$parm = i;df$label ='RQ';return(df) } ))
  y_sim_SE = do.call('rbind',lapply(1:nl, function(i){df = data.frame(t = t,y = mvrnorm(1, rep(0, n), Sigma_SE[[i]] )); df$parm = i;df$label ='SE';return(df) } ))
  y_sim_PR = do.call('rbind',lapply(1:nl, function(i){df = data.frame(t = t,y = mvrnorm(1, rep(0, n), Sigma_PR[[i]] )); df$parm = i;df$label ='PR';return(df) } ))
  y_sim_locPR = do.call('rbind',lapply(1:nl, function(i){df = data.frame(t = t,y = mvrnorm(1, rep(0, n), Sigma_locPR[[i]] )); df$parm = i;df$label ='locPR';return(df) } ))
  toplot_y_st_list[[s]] =  do.call('rbind',list(y_sim_SE,y_sim_RQ,y_sim_PR,y_sim_locPR))
  toplot_y_st_list[[s]]$s = s
}
toplot_y_st = do.call('rbind',toplot_y_st_list)

gg = ggplot(toplot_y_st[toplot_y_st$s %in% 1:3,], aes(x = t, y = y,color = factor(s)))
gg = gg + facet_grid(parm~label, scales = 'fixed' )
gg = gg + geom_line()
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'none')
gg = gg + labs(x = NULL, y = NULL)
ggsave(paste(dir,'/y_sim_stationary.png',sep = ''),width =10,height = 8)
ggsave(paste(dir,'/y_sim_stationary.eps',sep = ''),width =10,height = 8)


#### CONTOURS
ind_list = data.frame(t1 = t[1], t2 = t[seq(2,n,by = 2)],ind2 = seq(2,n,by = 2) )
ncases = nrow(ind_list)
params_list = unique(toplot_y_st[,c("parm" , "label")])

x = seq(-3,3,0.1)
XY = expand.grid(x,x)

df_SE  = apply( params_list[params_list $label =='SE',],1 ,function(rr){ do.call('rbind',
                                                                                 lapply(1:ncases, function(ss){ df = data.frame( label ='SE', parm = rr[1],  Y1 = XY[,1] ,Y2 = XY[,2], 
                                                                                                                                prob = dmvnorm(XY, mean = rep(0,2), sigma = Sigma_SE[[as.numeric(rr[1])]][c(1,ind_list[ss,3]),c(1,ind_list[ss,3])] ),
                                                                                                                                point_label = paste('t=',format(ind_list[ss,1:2], digits = 2), sep='', collapse=', ' )     )}    )   )}  )
df_SE = do.call('rbind',df_SE)

l_labels = paste('l=',l_list)
l_alpha_labels = apply(l_alpha_list,1,function(rr) paste(c('l=','a='),format(rr, digits = 2), sep='', collapse=', ' ))

df_SE$parm = factor(df_SE$parm)
levels(df_SE$parm) = l_labels 
gg = ggplot(data=df_SE,aes(x=Y1,y=Y2,z=prob,colour = parm)) 
gg = gg + facet_grid(point_label~label + parm, scales = 'fixed' )
gg = gg + geom_contour()
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'none')
ggsave(paste(dir,'/y_cont_stationary_SE.png',sep = ''),width =10,height = 8)
ggsave(paste(dir,'/y_cont_stationary_SE.eps',sep = ''),width =10,height = 8)


df_RQ  = apply( params_list[params_list $label =='RQ',],1 ,function(rr){ do.call('rbind',
                                                                                 lapply(1:ncases, function(ss){ df = data.frame( label ='RQ', parm = rr[1],  Y1 = XY[,1] ,Y2 = XY[,2], 
                                                                                                                                 prob = dmvnorm(XY, mean = rep(0,2), sigma = Sigma_RQ[[as.numeric(rr[1])]][c(1,ind_list[ss,3]),c(1,ind_list[ss,3])] ),
                                                                                                                                 point_label = paste('t=',format(ind_list[ss,1:2], digits = 2), sep='', collapse=', ' )     )}    )   )}  )
df_RQ = do.call('rbind',df_RQ)

df_RQ$parm = factor(df_RQ$parm)
levels(df_RQ$parm) = l_alpha_labels
gg = ggplot(data=df_RQ,aes(x=Y1,y=Y2,z=prob,colour = parm)) 
gg = gg + facet_grid(point_label~label + parm, scales = 'fixed' )
gg = gg + geom_contour()
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'none')
ggsave(paste(dir,'/y_cont_stationary_RQ.png',sep = ''),width =10,height = 8)
ggsave(paste(dir,'/y_cont_stationary_RQ.eps',sep = ''),width =10,height = 8)

l_alpha_labels = apply(l_alpha_list,1,function(rr) paste(c('l=','a='),format(rr*c(1,1/5), digits = 2), sep='', collapse=', ' ))

df_PR  = apply( params_list[params_list $label =='PR',],1 ,function(rr){ do.call('rbind',
                                                                                 lapply(1:ncases, function(ss){ df = data.frame( label ='PR', parm = rr[1],  Y1 = XY[,1] ,Y2 = XY[,2], 
                                                                                                                                 prob = dmvnorm(XY, mean = rep(0,2), sigma = Sigma_PR[[as.numeric(rr[1])]][c(1,ind_list[ss,3]),c(1,ind_list[ss,3])] ),
                                                                                                                                 point_label = paste('t=',format(ind_list[ss,1:2], digits = 2), sep='', collapse=', ' )     )}    )   )}  )
df_PR = do.call('rbind',df_PR)
df_PR$parm = factor(df_PR$parm)
levels(df_PR$parm) = l_alpha_labels
gg = ggplot(data=df_PR,aes(x=Y1,y=Y2,z=prob,colour = parm)) 
gg = gg + facet_grid(point_label~label + parm, scales = 'fixed' )
gg = gg + geom_contour()
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'none')
ggsave(paste(dir,'/y_cont_stationary_PR.png',sep = ''),width =10,height = 8)
ggsave(paste(dir,'/y_cont_stationary_PR.eps',sep = ''),width =10,height = 8)


df_locPR  = apply( params_list[params_list $label =='locPR',],1 ,function(rr){ do.call('rbind',
                                                                                 lapply(1:ncases, function(ss){ df = data.frame( label ='locPR', parm = rr[1],  Y1 = XY[,1] ,Y2 = XY[,2], 
                                                                                                                                 prob = dmvnorm(XY, mean = rep(0,2), sigma = Sigma_locPR[[as.numeric(rr[1])]][c(1,ind_list[ss,3]),c(1,ind_list[ss,3])] ),
                                                                                                                                 point_label = paste('t=',format(ind_list[ss,1:2], digits = 2), sep='', collapse=', ' )     )}    )   )}  )
df_locPR = do.call('rbind',df_locPR)
df_locPR$parm = factor(df_locPR$parm)
levels(df_locPR$parm) = l_alpha_labels
gg = ggplot(data=df_locPR,aes(x=Y1,y=Y2,z=prob,colour = parm)) 
gg = gg + facet_grid(point_label~label + parm, scales = 'fixed' )
gg = gg + geom_contour()
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'none')
ggsave(paste(dir,'/y_cont_stationary_locPR.png',sep = ''),width =10,height = 8)
ggsave(paste(dir,'/y_cont_stationary_locPR.eps',sep = ''),width =10,height = 8)



######### PREDICTION
id_obs = sample(1:n,4,replace= F)
nobs = length(id_obs)
npred = n - nobs
t_obs = t[id_obs]
t_pred = t[-id_obs]
y_obs = subset(toplot_y_st, toplot_y_st$s == 1 & toplot_y_st$label=='RQ' & toplot_y_st$parm == 1 & toplot_y_st$t %in% t_obs)[,'y']

mu_pred_SE = do.call('rbind',lapply(1:nl, function(i){df = data.frame(t = c(t_obs,t_pred),
                                                                      y = c(y_obs,Sigma_SE[[i]][-id_obs,id_obs] %*% solve(Sigma_SE[[i]][id_obs,id_obs]) %*% matrix(y_obs,ncol = 1)),
                                                                      conf_spread = c(rep(0,nobs), diag(Sigma_SE[[i]][-id_obs,-id_obs] - Sigma_SE[[i]][-id_obs,id_obs] %*% solve( Sigma_SE[[i]][id_obs,id_obs]) %*% Sigma_SE[[i]][id_obs,-id_obs]) ),
                                                                      isObs = c(rep(1,nobs),rep(0,npred) ),    
                                                                      parm =  i, 
                                                                      label ='SE');return(df) } ))

mu_pred_RQ = do.call('rbind',lapply(1:nla, function(i){df = data.frame(t = c(t_obs,t_pred),
                                                                      y = c(y_obs,Sigma_RQ[[i]][-id_obs,id_obs] %*% solve(Sigma_RQ[[i]][id_obs,id_obs]) %*% matrix(y_obs,ncol = 1)),
                                                                      conf_spread = c(rep(0,nobs), diag(Sigma_RQ[[i]][-id_obs,-id_obs] - Sigma_RQ[[i]][-id_obs,id_obs] %*% solve( Sigma_RQ[[i]][id_obs,id_obs]) %*% Sigma_RQ[[i]][id_obs,-id_obs]) ),
                                                                      isObs = c(rep(1,nobs),rep(0,npred) ),    
                                                                      parm =  i, 
                                                                      label ='RQ');return(df) } ))



mu_pred_PR = do.call('rbind',lapply(1:nl, function(i){df = data.frame(t = c(t_obs,t_pred),
                                                                      y = c(y_obs,Sigma_PR[[i]][-id_obs,id_obs] %*% solve(Sigma_PR[[i]][id_obs,id_obs]) %*% matrix(y_obs,ncol = 1)),
                                                                      conf_spread = c(rep(0,nobs), diag(Sigma_PR[[i]][-id_obs,-id_obs] - Sigma_PR[[i]][-id_obs,id_obs] %*% solve( Sigma_PR[[i]][id_obs,id_obs]) %*% Sigma_PR[[i]][id_obs,-id_obs]) ),
                                                                      isObs = c(rep(1,nobs),rep(0,npred) ),    
                                                                      parm =  i, 
                                                                      label ='PR');return(df) } ))

mu_pred_locPR = do.call('rbind',lapply(1:nl, function(i){df = data.frame(t = c(t_obs,t_pred),
                                                                      y = c(y_obs,Sigma_locPR[[i]][-id_obs,id_obs] %*% solve(Sigma_locPR[[i]][id_obs,id_obs]) %*% matrix(y_obs,ncol = 1)),
                                                                      conf_spread = c(rep(0,nobs), diag(Sigma_locPR[[i]][-id_obs,-id_obs] - Sigma_locPR[[i]][-id_obs,id_obs] %*% solve( Sigma_locPR[[i]][id_obs,id_obs]) %*% Sigma_locPR[[i]][id_obs,-id_obs]) ),
                                                                      isObs = c(rep(1,nobs),rep(0,npred) ),    
                                                                      parm =  i, 
                                                                      label ='locPR');return(df) } ))

toplot_y_pred_st =  do.call('rbind',list(mu_pred_SE,mu_pred_RQ,mu_pred_PR,mu_pred_locPR))
toplot_y_pred_st$conf_up = toplot_y_pred_st$y + 1.96 * toplot_y_pred_st$conf_spread
toplot_y_pred_st$conf_down = toplot_y_pred_st$y - 1.96 * toplot_y_pred_st$conf_spread
gg = ggplot(data = toplot_y_pred_st)
gg = gg + facet_grid(parm~label, scales = 'fixed' )
gg = gg + geom_line(aes(x = t, y = y ) )
gg = gg + geom_ribbon(aes(x = t, ymin = conf_down,ymax = conf_up),alpha = 0.5)
gg = gg + geom_point(data = data.frame(y_obs = y_obs, t_obs = t_obs), aes(y = y_obs, x = t_obs) ,color ='red', size = 3)
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
#gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'none')
gg = gg + labs(x = NULL, y = NULL)
ggsave(paste(dir,'/y_pred_stationary.png',sep = ''),width =10,height = 8)
ggsave(paste(dir,'/y_pred_stationary.eps',sep = ''),width =10,height = 8)



confidence_interval = sigma * 1.96
q_80 = qmvnorm(p = 0.95, sigma = function_Sigma_SE(1)[1:2,1:2])



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


