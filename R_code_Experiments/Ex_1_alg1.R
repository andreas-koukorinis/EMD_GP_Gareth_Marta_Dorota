library(mvtnorm)
library(fBasics)
library(RConics)   #optional - for plot
library(reshape2)  #optional - for plot
library(ggplot2)   #optional - for plot


set.seed(1)

i= 1
j= 1
N= 100
m= 1
kernel_family= 0 


if(kernel_family == 0){
  
  Omega_true = c(0.25,0.5, 3,4,5, 10, 50)
  
  Omega_test = c(seq(0.01,10,0.01), 15,20,30,40,50,60)
  Omega_test = unique(c( Omega_true, Omega_test))    
  
  psi_test = Omega_test[i]
  psi_true = Omega_true[j]
  
  k = function(x1,x2,l ) {exp(- (x1-x2)^2 / (2 * l^2) ) } 
  
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
  
} else if(kernel_family == 2){        
  
  l_list_true= c(0.25,0.5, 3,4,5, 10, 50)
  alpha_list_true= c(0.25, 0.5, 1,2)
  Omega_true = expand.grid(l_list_true,alpha_list_true)
  
  l_list_test= c(seq(0.01,10,0.01), 15,20,30,40,50,60)
  alpha_list_test= c(0.25, 0.5, 1,2)
  l_list_test= unique(c(l_list_true, l_list_test))
  alpha_list_test= unique(c(alpha_list_true, alpha_list_test))
  Omega_test= expand.grid(l_list_test, alpha_list_test)
  
  Omega_test = unique(rbind(Omega_true, Omega_test), MARGIN = 1)   
  
  psi_test = Omega_test[i,]
  psi_true = Omega_true[j,]
  
  k = function(x1,x2,l ,alpha ) { ( 1 + abs(x1-x2)^2/(2 * alpha * l^2)  )^{-alpha} } #NO ABS RIGHT?
  
} else if(kernel_family == 3){       
  
  l_list_true= c(0.25,0.5, 3,4,5, 10, 50)
  p_list_true= c(0.25, 0.5, 1,2)
  Omega_true = expand.grid(l_list_true,p_list_true)
  
  l_list_test= c(seq(0.01,10,0.01), 15,20,30,40,50,60)
  p_list_test= c(0.25, 0.5, 1,2)
  l_list_test= unique(c(l_list_true, l_list_test))
  p_list_test= unique(c(p_list_true, p_list_test))
  Omega_test= expand.grid(l_list_test, p_list_test)
  
  Omega_test = unique(rbind(Omega_true, Omega_test), MARGIN = 1)   
  
  psi_test = Omega_test[i,]
  psi_true = Omega_true[j,]
  
  k = function(x1,x2,l ,p) {{exp(-2 * sin(pi*(x1 - x2) / p)^2 / (l^2) ) }  * {exp(- (x1-x2)^2 / (2 * l^2) ) } } 
  
} else if(kernel_family == 4){        
  
  Omega_true = c(0.25,0.5, 3,4,5, 10, 50)   
  
  Omega_test = c(seq(0.01,10,0.01), 15,20,30,40,50,60)
  Omega_test = unique(c( Omega_true, Omega_test))  
  
  psi_test = Omega_test[i]
  psi_true = Omega_true[j]
  
  k = function(x1,x2,sigma2 ){x1 * x2 + sigma2} 
  
} else if(kernel_family == 5){        #MATERN3 or MATERN5 - PLOT THEM  DOROTA!!!!! 
  
  sigmaf_list_true= c(0.25,0.5, 3,4,5, 10, 50)
  sigmal_list_true= c(0.25,0.5, 3,4,5, 10, 50)
  r_list_true= c(0.25,0.5, 3,4,5, 10, 50)
  Omega_true = expand.grid(sigmaf_list_true,sigmal_list_true, r_list_true)
  
  sigmaf_list_test= c(0.25,0.5, 3,4,5, 10, 50)
  sigmal_list_test= c(0.25,0.5, 3,4,5, 10, 50)
  r_list_test= c(0.25,0.5, 3,4,5, 10, 50)
  sigmaf_list_test= unique(c(sigmaf_list_true, sigmaf_list_test))
  sigmal_list_test= unique(c(sigmal_list_true, sigmal_list_test))
  r_list_test = unique(c(r_list_true, r_list_test))
  Omega_test= expand.grid(sigmaf_list_test,sigmal_list_test, r_list_test)
  
  Omega_test = unique(rbind(Omega_true, Omega_test), MARGIN = 1)   
  
  psi_test = Omega_test[i,]
  psi_true = Omega_true[j,]
  
  k = function(x1,x2,sigmaf, sigmal, r) {sigmaf*(1 + ((sqrt(3)*r)/ sigmal ) ) * exp(-(sqrt(3)*r)/sigmal)} 
  
}



t= seq(0,1, length.out = N)

K_true= matrix(0, N,N) 
K_test= matrix(0, N,N)


for(n1 in (1:N)){
  
  for(n2 in (n1:N)){
    
    K_true[n1,n2]= k(t[n1],t[n2], l = psi_true )  
    K_test[n1,n2]= k(t[n1],t[n2], l = psi_test)
    
  }
}


K_true = K_true + t(K_true)
diag(K_true) = diag(K_true) / 2


K_test = K_test + t(K_test)
diag(K_test) = diag(K_test) / 2


y_m= rmvnorm(1,rep(0,N), K_true)
y_m= y_m /as.numeric(sqrt(var(t(y_m))))

S_m= t(y_m)%*%y_m

a_ijm= 0

a_ijm= 2 - norm((S_m/norm( S_m , type = "F") -  K_test/norm( K_test , type = "F")) , type = "F")

return_list = list(a_ijm = a_ijm, i = i, j = j, m = m , N = N, psi_true = psi_true, psi_test = psi_test)



#DOROTA
#TRY the trace for CKTA - check which operation is optimal CHECK

toplot = melt(adjoint(K_test))
gg = ggplot(toplot, aes(x = Var1, y = Var2,fill = value))
gg = gg + geom_tile()
gg = gg + scale_fill_gradientn(colors = c('lightyellow','blue','red'))
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'bottom')
gg = gg + labs(x = NULL, y = NULL)








