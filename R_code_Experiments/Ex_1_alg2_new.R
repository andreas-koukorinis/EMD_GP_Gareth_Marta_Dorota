library(mvtnorm)
library(ggplot2)
library(dplyr)

set.seed(1)

N= c(10^2,10^3)
M = 10
t<- lapply(N, function(n) seq(0,1, length.out = n) )

kernel_family=0 

if(kernel_family == 0){
  
  k =function(x1,x2,l ) {exp(- (x1-x2)^2 / (2 * l^2) ) }
  #KSE
  Omega_true = c(0.25,0.5, 3,4,5, 10, 50)
  
  Omega_test = c(seq(0.01,1,0.1), 15,20,30,40,50,60)
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

####### DEFINING K MATRICES SKELETONS
K_true =  lapply(N, function(n) lapply(1:J, function(s) matrix(0, n, n) ) )### kernel func evaluated on deterministic vector of times 't' - does not depend on 'm'
K_test =lapply(N, function(n) lapply(1:I, function(i) matrix(0,n,n) ))  ### kernel func evaluated on deterministic vector of times 't' - does not depend on 'm'

for( n in 1:length(t)){
  for (j in 1:J) {
    psi_true = Omega_true[j]
      for(n1 in 1:N[n] ){
        for(n2 in n1: N[n]){
          K_true[[n]][[j]][n1,n2]=k(t[[n]][n1],t[[n]][n2], l = psi_true )  
      }
    }
    K_true[[n]][[j]] = K_true[[n]][[j]] + t(K_true[[n]][[j]])
    diag(K_true[[n]][[j]]) = diag(K_true[[n]][[j]])/2
  }
}

for( n in 1:length(t)){
  for (j in 1:I) {
    psi_test = Omega_test[j]
    for(n1 in 1:N[n]){
      for(n2 in n1:N[n]){
        K_test[[n]][[j]][n1,n2]=k(t[[n]][n1],t[[n]][n2], l = psi_test )
      }
    }
    K_test[[n]][[j]]=K_test[[n]][[j]] + t(K_test[[n]][[j]])
    diag(K_test[[n]][[j]]) =diag(K_test[[n]][[j]])/2   
  }
}


#### Simulating Gaussian Processes with true Gram Matrices K_true
y_m = list() ### correct, the random vector y is depndent on m since is mth simulation from a Gaussian distribution with the cov matrix == gram matrix K_true
for( n in 1:length(t)){
  y_m[[n]] = vector(mode = "list", length = J)
  for(j in (1:J)){
    for(m in (1:M)){
      y_m[[n]][[j]][[m]] = rmvnorm(1,rep(0,N[n]),  K_true[[n]][[j]])
      y_m[[n]][[j]][[m]] = y_m[[n]][[j]][[m]] / sd(y_m[[n]][[j]][[m]]) ## no need for an additional package for 'sd' fun
    }
  }
}

S_m = list()
for( n in 1:length(t)){
  S_m[[n]]= vector(mode = "list", length = J)
  for(j in (1:J)){
    for(m in (1:M))
      S_m[[n]][[j]][[m]]=t(y_m[[n]][[j]][[m]])%*%y_m[[n]][[j]][[m]]
  }
}


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

calculate_CKTA_old  = function(Kx,Ky){
  F_norm_Kx = Kx/norm(Kx, type = "F")
  F_norm_Ky = Ky/norm(Ky, type = "F")
  ckta =  2 - norm( F_norm_Kx  - F_norm_Ky , type = "F")
  return(ckta)
}


##### TO ANALYSE FOR MARTA
# j = 1
# i = 1
# m = 1
# calculate_CKTA_old(Kx =S_m[[j]][[m]] ,Ky = K_test[[i]] )
# calculate_CKTA(Kx = S_m[[j]][[m]] ,Ky = K_test[[i]] )
# 
# calculate_CKTA_old(Kx = K_test[[i]] ,Ky =  K_test[[i]] )
# calculate_CKTA(Kx =  K_test[[i]] ,Ky = K_test[[i]] )
# 
# calculate_CKTA_old(Kx = -K_test[[i]] ,Ky =  K_test[[i]] )
# calculate_CKTA(Kx = -K_test[[i]] ,Ky =  K_test[[i]] )
# 
# calculate_CKTA_old(Kx =diag(N) ,Ky = diag(N) )
# calculate_CKTA(Kx = diag(N) ,Ky = diag(N) )


a_ijm = t( apply(CKTA_grid, 1, function(irow){ n = irow[1]; j = irow[2]; i = irow[3]; m = irow[4];  c( calculate_CKTA_old( Kx = S_m[[n]][[j]][[m]] , Ky = K_test[[n]][[i]] ) , n,j,i,m,  Omega_true[j],  Omega_test[i])  } ) )
colnames(a_ijm) = c('a_ijm','N','j','i','m','psi_true','psi_test')

a_ijm_true_vs_true = a_ijm[ a_ijm[,3] == a_ijm[,4], ]
colnames(a_ijm_true_vs_true)[1] = 'a_jjm'
a_ijm_final = merge(a_ijm, a_ijm_true_vs_true[,-c(4,7)], all.x = TRUE)



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
gg = gg + geom_vline( aes(xintercept = psi_true), linetype = 2 )
gg = gg + geom_boxplot(alpha = 0, fill = NA, color = 'red')
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + labs(title = 'log of square diff between CKTA with true Psi and best Psi')
gg = gg + theme(legend.position = 'bottom',legend.text = element_text(size = 10))
gg = gg + theme(axis.text.x = element_text(angle = 90))
gg = gg + guides(col = guide_legend(nrow = 2))
gg


gg = ggplot(a_ijm_optimal, aes(x = psi_test))
gg = gg + facet_grid( N_l ~ psi_true_l,  scales = 'free_y')
gg = gg + geom_bar()
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + labs(title = 'optimal Psi test choice - counts')
gg = gg + theme(legend.position = 'bottom',legend.text = element_text(size = 10))
gg = gg + theme(axis.text.x = element_text(angle = 90))
gg = gg + guides(col = guide_legend(nrow = 2))
gg








