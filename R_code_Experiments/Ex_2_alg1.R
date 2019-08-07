library("EMD")
library("splines")
library("seewave")
library("hht")
library(mvtnorm)
library("phonTools")

set.seed(1)

m=1
N=1000   #1000
kernel_family= 1 


if(kernel_family == 0){
  Omega = c(0.05, 0.10, 0.15, 0.25,0.5, 1.5,2,2.5, 3, 5, 10,30)
 
  k = function(x1,x2,l ) {exp(-(x1-x2)^2 / (2 * l^2) ) }
  
} else if(kernel_family == 1){
  l_list= c(0.05, 0.10, 0.15, 0.25, 0.5, 1)
  p_list= c(0.25, 0.5, 1)
  Omega = expand.grid(l_list,p_list)
  
  k= function(x1,x2,l ,p ){exp(-2 * sin(pi*(x1 - x2) / p)^2 / (l^2) ) }
  
} 

J= dim(Omega)[1]    #3
#j= dim(Omega)[2]
t= seq(0,1, length.out = N)

K= lapply(1:J, function(h)matrix(0,N ,N)) 


for (h in 1:J) {
  for(n1 in (1:N)){
    for(n2 in (n1:N)){
      
      psi = Omega[h,]
      K[[h]][n1,n2]= k(t[n1],t[n2], l = psi[1], p = psi[2] )    
     
    }
  }
}


for (h in 1:J) {

  K[[h]]= K[[h]] + t(K[[h]])
  diag(K[[h]]) = diag(K[[h]])/2
  
}



y_m= vector(mode = "list", length = J)

for(h in (1:J)){

    y_m[[h]]= rmvnorm(1,rep(0,(dim(K[[h]])[1])),  K[[h]])

}


#Sum over all the hyperparameters
x_m= Reduce('+', y_m)

#Fit a spline
x_m_spl= splinefun(t, x_m)


hat_x_m= x_m_spl(t, deriv = 0)   #knots same as obs. points

#EMD
dec= emd(hat_x_m)


IF= vector(mode = "list", length = dim(dec$imf)[2])

for (i in 1: ( dim(dec$imf)[2])) {
  IF[[i]]= ifreq(dec$imf[,i], f = 4000)$f[,2]
  
}


#SPECTROGRAM
spectrogram(hat_x_m, fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000)

spectrogram(dec$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000)

spectrogram(dec$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000)

spectrogram(dec$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -1000, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000)

spectrogram(dec$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000)

spectrogram(dec$imf[,5], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000)




