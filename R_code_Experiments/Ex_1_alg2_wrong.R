library(mvtnorm)
library(fBasics)

set.seed(1)

N=1000
m=N/100

tt=seq(0,1, length.out = N) # USE M AS SAME NOTATION IN OUR LATEX FILE - M AND NOT H
t=vector(mode = "list", length = m)      # TAKE SUBSET OF THE GRAM MATRIX  

for (i in 1:m) {
  t[[i]]=tt[1:(m*10*i)]
}


kernel_family=0 

#KSE
Omega_true = c(0.25,0.5, 3,4,5, 10, 50)

Omega_test = c(seq(0.01,10,0.01), 15,20,30,40,50,60)
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



K_true=lapply(1:J, function(j)lapply(1:m,function(i)matrix(0,(m*10*i) ,(m*10*i))) )
K_test=lapply(1:I, function(j)lapply(1:m,function(i)matrix(0,(m*10*i) ,(m*10*i))) )


# K_true=lapply(1:J, function(j)lapply(1:m,function(i)matrix(0,(100) ,(100))) )
# K_test=lapply(1:I, function(j)lapply(1:m,function(i)matrix(0,(100) ,(100))) )


for (j in 1:J) {
  for(h in 1:m){
    for(n1 in (1:(dim(K_true[[j]][[h]])[1]))){
      
      for(n2 in (n1:(dim(K_true[[j]][[h]])[1]))){
        
        psi_true = Omega_true[j]
        K_true[[j]][[h]][n1,n2]=k(t[[h]][n1],t[[h]][n2], l = psi_true )  
        
      }
    }
  }
}


for (j in 1:I) {
  for(h in 1:m){
    for(n1 in (1:(dim(K_test[[j]][[h]])[1]))){
      
      for(n2 in (n1:(dim(K_test[[j]][[h]])[1]))){
        
        psi_test = Omega_test[j]
        K_test[[j]][[h]][n1,n2]=k(t[[h]][n1],t[[h]][n2], l = psi_test )  
        
      }
    }
  }
}

#DO EVERYTHING IN ONE LOOP ABOVE
for (j in 1:J) {
  for(h in 1:m){
    
        K_true[[j]][[h]]=K_true[[j]][[h]] + t(K_true[[j]][[h]])
        diag(K_true[[j]][[h]]) =diag(K_true[[j]][[h]])/2

  }
}


for (j in 1:I) {
  for(h in 1:m){
    
      K_test[[j]][[h]]=K_test[[j]][[h]] + t(K_test[[j]][[h]])
      diag(K_test[[j]][[h]]) =diag(K_test[[j]][[h]])/2        

  }
}



y_m=vector(mode = "list", length = J)

for(j in (1:J)){
for(h in (1:m)){
  
  y_m[[j]][[h]]=rmvnorm(1,rep(0,(dim(K_true[[j]][[h]])[1])),  K_true[[j]][[h]])
  y_m[[j]][[h]]=y_m[[j]][[h]] /as.numeric(sqrt(var(t(y_m[[j]][[h]]))))
}
}


S_m=lapply(1:J, function(j)lapply(1:m,function(i)matrix(0,(m*10*i) ,(m*10*i))) )
#S_m=lapply(1:J, function(j)lapply(1:m,function(i)matrix(0,(100) ,(100))) )

for(j in (1:J)){
  for(h in (1:m))
    
    S_m[[j]][[h]]=t(y_m[[j]][[h]])%*%y_m[[j]][[h]]
  
}



F_norm_S_m=vector(mode = "list", length = J)
for(j in (1:J)){
  for(h in (1:m))
    
    F_norm_S_m[[j]][[h]]=(S_m[[j]][[m]]/norm(S_m[[j]][[m]], type = "F"))
  
}

F_norm_K_test=vector(mode= "list", length = I)
for(i in (1:I)){
  for(h in (1:m))
    
    F_norm_K_test[[i]][[h]]=(K_test[[i]][[m]]/norm(K_test[[i]][[m]], type = "F"))
  
}


a_ijm=lapply(1:J, function(j){
  
  a_ijm1=lapply(1:I, function(i){
    
    a_ijm2=sapply(1:m, function(h){
      
      2 - norm( F_norm_S_m[[j]][[h]] - F_norm_K_test[[i]][[h]] , type = "F")
      
    })
    
  })
  
  
})



return_list=lapply(1:J, function(j){
  
  return_list1=lapply(1:I, function(i){
    
    return_list2=sapply(1:m, function(h){
      
      c(a_ijm = a_ijm[[j]][[i]][h], i = i, j = j, m = h ,  
           psi_true = Omega_true[j], psi_test = Omega_test[i]) #N = N,
      
    })
    
  })
  
  
})



#return_list = list(a_ijm = a_ijm[[j]][[i]][[h]], i = i, j = j,
#                   m = m , N = N, psi_true = psi_true, psi_test = psi_test)





#DOROTA - BOX-PLOT EXERCISE









