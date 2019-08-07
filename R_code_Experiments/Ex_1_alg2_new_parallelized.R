library(mvtnorm)
library(foreach)
library(doParallel)

dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)
set.seed(1)
n.cores = detectCores() - 1 ## to ensure one core for the operating system not to crash

N=c(10^2,10^3)
M = 10^2
t = lapply(N, function(n) seq(0,1, length.out = n))

kernel_family=0 
if(kernel_family == 0){
  
  k =function(x1,x2,l ) {exp(- (x1-x2)^2 / (2 * l^2) ) }
  #KSE
  Omega_true = c(0.05, 0.25, 0.5, 3,5,10)
  
  Omega_test = c(seq(0.01,1,0.1), 15,20,30)
  Omega_test = unique(c( Omega_true, Omega_test))    
  
  J=length(Omega_true) #M_true
  I=length(Omega_test) #M_test
  
} else if(kernel_family == 1){         
  
  l_list_true= c(0.25,0.5, 3,4,5, 10, 50)
  p_list_true= c(0.25, 0.5, 1,2)
  Omega_true = expand.grid(l_list_true,p_list_true)
  
  l_list_test= c(seq(0.01,10,0.01), 15,20,30,40,50,60)
  p_list_test= c(0.25, 0.5, 1,2)
  l_list_test = unique(c(l_list_true, l_list_test))
  p_list_test = unique(c(p_list_true, p_list_test))
  Omega_test= expand.grid(l_list_test, p_list_test)
  
  Omega_test = unique(rbind(Omega_true, Omega_test), MARGIN = 1)   
  
  psi_test = Omega_test[i,]
  psi_true = Omega_true[j,]
  
  k= function(x1,x2,l ,p ){exp(-2 * sin(pi*(x1 - x2) / p)^2 / (l^2) ) } 
  
}


calculate_gram_mat = function(N0,t0,k_fun){
  K = matrix(0, N0, N0)
  for(n1 in 1:N0 ){
    for(n2 in n1: N0){
      K[n1,n2]=k_fun(t0[n1],t0[n2])  
    }
  }
  K = K + t(K)
  diag(K) = diag(K)/2
  return(K)
}

###### CALCULATING GRAM MATRICES
cl = makeForkCluster(n.cores) ### THIS ONE WORKS ON IOS AND LINUX ONLY< FOR MARTA TO FIGURE OUT WHICH ONE WOULD WORK ON WINDOWS
registerDoParallel(cl) 
K_true = foreach(n = 1:length(N)) %:%
  foreach(j = 1:J) %dopar% {
    K = calculate_gram_mat(N0 = N[n],t0 = t[[n]],k_fun = function(t1,t2) k(t1,t2, l = Omega_true[j] ) )
  }
stopCluster(cl)

cl = makeForkCluster(n.cores)
registerDoParallel(cl) 
K_test = foreach(n = 1:length(N)) %:%
  foreach(i = 1:I) %dopar% {
    K = calculate_gram_mat(N0 = N[n],t0 = t[[n]],k_fun = function(t1,t2) k(t1,t2, l = Omega_test[i] ) )
  }
stopCluster(cl)


######## CALCULATING CENTRAL KERNEL TARGET ALIGMENT
CKTA_grid = expand.grid(1:length(N), 1:J,1:I,1:M)

calculate_CKTA = function(Kx,Ky){
  Kxc = t(t(Kx)-colMeans(Kx))
  Kxc = Kxc-rowMeans(Kxc)
  
  Kyc = t(t(Ky)-colMeans(Ky))
  Kyc = Kyc-rowMeans(Kyc)
  
  numer <- sum(Kx*Kyc)
  denomx <- sum(Kx*Kxc)
  denomy <- sum(Ky*Kyc)
  return(  numer/sqrt(denomx*denomy) )
}


compute_parallel_CKTA = function(ii){
  irow = as.numeric(CKTA_grid[ii,])
  n = irow[1]
  j = irow[2]
  i = irow[3]
  m = irow[4]
  y0 = rmvnorm(1,rep(0,N[n]),  K_true[[n]][[j]])
  #y0 = y0 / sd(y0)
  return( c( calculate_CKTA( Kx = t(y0)%*%y0 , Ky = K_test[[n]][[i]] ) , n,j,i,m,  Omega_true[j],  Omega_test[i]) )
}

library(parallel)  
a_ijm_list = mclapply(1:nrow(CKTA_grid), compute_parallel_CKTA, mc.cores = 9)

a_ijm = do.call('rbind',a_ijm_list)
colnames(a_ijm) = c('a_ijm','N','j','i','m','psi_true','psi_test')

a_ijm_true_vs_true = a_ijm[ a_ijm[,3] == a_ijm[,4], ]
colnames(a_ijm_true_vs_true)[1] = 'a_jjm'
a_ijm_final = merge(a_ijm, a_ijm_true_vs_true[,-c(4,7)], all.x = TRUE)

library(ggplot2)
library(dplyr)

a_ijm_optimal = a_ijm_final %>% group_by(N,j,m) %>% top_n(n = 1, wt = a_ijm)  %>% arrange(j)
a_ijm_optimal$diff_ckta = log((a_ijm_optimal$a_ijm - a_ijm_optimal$a_jjm)**2)
a_ijm_optimal$psi_test = factor(a_ijm_optimal$psi_test, levels = sort(Omega_test))
a_ijm_optimal$psi_true_l = factor(a_ijm_optimal$psi_true, levels = Omega_true)
levels(a_ijm_optimal$psi_true_l) =  paste( 'psi_true:', levels(a_ijm_optimal$psi_true_l))
a_ijm_optimal$N = factor(a_ijm_optimal$N, levels = 1:length(N))
a_ijm_optimal$N_l = a_ijm_optimal$N
levels(a_ijm_optimal$N_l) = paste('N:', N)

gg = ggplot(a_ijm_optimal, aes(x = psi_test,y = diff_ckta))
gg = gg + facet_grid( N_l ~ psi_true_l,  scales = 'free_y')
gg = gg + geom_jitter(aes(color = psi_test),alpha = 0.5,width = 0.05,shape = 16) 
gg = gg + geom_boxplot(alpha = 0, fill = NA, color = 'red')
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + labs(title = 'log of square diff between CKTA with true Psi and best Psi, nsim = 100')
gg = gg + theme(legend.position = 'bottom',legend.text = element_text(size = 10))
gg = gg + theme(axis.text.x = element_text(angle = 90))
gg = gg + guides(col = guide_legend(nrow = 2))
ggsave(paste(dir,'/CKTA_SE_no_sd.png',sep = ''), width =10,height = 7)



gg = ggplot(a_ijm_optimal, aes(x = psi_test))
gg = gg + facet_grid( N_l ~ psi_true_l,  scales = 'free_y')
gg = gg + geom_bar()
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + labs(title = 'Optimal Psi test choice - counts, nsim = 100')
gg = gg + theme(legend.position = 'bottom',legend.text = element_text(size = 10))
gg = gg + theme(axis.text.x = element_text(angle = 90))
gg = gg + guides(col = guide_legend(nrow = 2))
ggsave(paste(dir,'/count_SE_no_sd.png',sep = ''), width =10,height = 7)





