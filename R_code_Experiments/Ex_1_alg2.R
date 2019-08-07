library(mvtnorm)
library(fBasics)

set.seed(1)

N=100
M = 10
t<- seq(0,1, length.out = N)

kernel_family=0 

#KSE
Omega_true = c(0.25,0.5, 3,4,5, 10, 50)

Omega_test = c(seq(0.01,1,0.1), 15,20,30,40,50,60)
Omega_test = unique(c( Omega_true, Omega_test))    

J=length(Omega_true) #M_true
I=length(Omega_test) #M_test

if(kernel_family == 0){
  
  k =function(x1,x2,l ) {exp(- (x1-x2)^2 / (2 * l^2) ) }
  
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


K_true=lapply(1:J, function(m)lapply(1:M,function(i)matrix(0,(N) ,(N))) )
K_test=lapply(1:I, function(m)lapply(1:M,function(i)matrix(0,(N) ,(N))) )


for (j in 1:J) {
  for(m in 1:M){
    for(n1 in (1:(dim(K_true[[j]][[m]])[1]))){
      
      for(n2 in (n1:(dim(K_true[[j]][[m]])[1]))){
        
        psi_true = Omega_true[j]
        K_true[[j]][[m]][n1,n2]=k(t[n1],t[n2], l = psi_true )  
       
      }
    }
  }
}


for (j in 1:I) {
  for(m in 1:M){
    for(n1 in (1:(dim(K_test[[j]][[m]])[1]))){
      
      for(n2 in (n1:(dim(K_test[[j]][[m]])[1]))){
        
        psi_test = Omega_test[j]
        K_test[[j]][[m]][n1,n2]=k(t[n1],t[n2], l = psi_test )
        
      }
    }
  }
}


for (j in 1:J) {
  for(m in 1:M){
    
    K_true[[j]][[m]] = K_true[[j]][[m]] + t(K_true[[j]][[m]])
    diag(K_true[[j]][[m]]) = diag(K_true[[j]][[m]])/2
    
  }
}


for (j in 1:I) {
  for(m in 1:M){
    
    K_test[[j]][[m]]=K_test[[j]][[m]] + t(K_test[[j]][[m]])
    diag(K_test[[j]][[m]]) =diag(K_test[[j]][[m]])/2        
    
  }
}



y_m=vector(mode = "list", length = J)

for(j in (1:J)){
  for(m in (1:M)){
    
    y_m[[j]][[m]]=rmvnorm(1,rep(0,(dim(K_true[[j]][[m]])[1])),  K_true[[j]][[m]])
    y_m[[j]][[m]]=y_m[[j]][[m]] /as.numeric(sqrt(var(t(y_m[[j]][[m]]))))
  }
}


S_m=lapply(1:J, function(m)lapply(1:M,function(i)matrix(0,(N) ,(N))) )

for(j in (1:J)){
  for(m in (1:M))
    
    S_m[[j]][[m]]=t(y_m[[j]][[m]])%*%y_m[[j]][[m]]
  
}



F_norm_S_m=vector(mode = "list", length = J)
for(j in (1:J)){
  for(m in (1:M))
    
    F_norm_S_m[[j]][[m]]=(S_m[[j]][[m]]/norm(S_m[[j]][[m]], type = "F"))
  
}

F_norm_K_test=vector(mode= "list", length = I)
for(i in (1:I)){
  for(m in (1:M))
    
    F_norm_K_test[[i]][[m]]=(K_test[[i]][[m]]/norm(K_test[[i]][[m]], type = "F"))
  
}


a_ijm=lapply(1:J, function(j){
  
  a_ijm1=lapply(1:I, function(i){
    
    a_ijm2=sapply(1:M, function(m){
      
      2 - norm( F_norm_S_m[[j]][[m]] - F_norm_K_test[[i]][[m]] , type = "F")
      
    })
    
  })
  
  
})



return_list=lapply(1:J, function(j){
  
  return_list1=lapply(1:I, function(i){
    
    return_list2=sapply(1:M, function(m){
      
      c(a_ijm = a_ijm[[j]][[i]][m], i = i, j = j, m = m ,  
        psi_true = Omega_true[j], psi_test = Omega_test[i], N = N) 
      
    })
    
  })
  
  
})





m_aijm = vector(mode = "list", length = J)
m_aijm1 = matrix(0, nrow = I, ncol = M)

difference_psi = vector(mode = "list", length = J)
difference_psi1 = matrix(0, nrow = I, ncol = M)

difference_ckta = vector(mode = "list", length = J)
difference_ckta1 = matrix(0, nrow = I, ncol = M)

psi_true1 = matrix(0, nrow = J, ncol = 1)

psi_test1 = matrix(0, nrow = I, ncol = 1)

for (j in 1:J) {
  for (i in 1:I) {
    for(m in 1:M){
      m_aijm1[i,m] =  return_list[[j]][[i]][1,m]
      
      dif_psi = return_list[[j]][[i]][6,1] - return_list[[j]][[i]][5,1]
      difference_psi1[i,m] =  dif_psi
      rm(dif_psi)
      
#      dif_ckta = 
      
    }
    
    psi_test1[i,1] = return_list[[1]][[i]][6,1]
  }
  
  
  m_aijm[[j]] = m_aijm1
  m_aijm1 = matrix(0, nrow = I, ncol = M)
  
  difference_psi[[j]] = difference_psi1
  difference_psi1 = matrix(0, nrow = I, ncol = M)
  
  psi_true1[j,1] = return_list[[j]][[1]][5,1]
  
}

rm(difference_psi1, m_aijm1)

