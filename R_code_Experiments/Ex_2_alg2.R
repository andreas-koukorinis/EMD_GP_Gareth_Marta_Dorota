library("EMD")
library("splines")
library("seewave")
library("hht")
library(mvtnorm)
library("phonTools")

set.seed(1)

m=10
N=100   #100
t= seq(0,1, length.out = N)

kernel_family= 1 

if(kernel_family == 0){
  Omega = c(0.05, 0.10, 0.15, 0.25,0.5, 1.5,2,2.5, 3, 5, 10,30)
  
  k = function(x1,x2,l ) {exp(-(x1-x2)^2 / (2 * l^2) ) }
  
} else if(kernel_family == 1){
  l_list= c(0.05, 0.10, 0.15, 0.2, 0.5, 1)
  p_list= c(0.25, 0.5, 1)
  Omega = expand.grid(l_list,p_list)
  
  k= function(x1,x2,l ,p ){exp(-2 * sin(pi*(x1 - x2) / p)^2 / (l^2) ) }
  
} 

Omega0 = rbind(Omega[1,],Omega[7,], Omega[13,])
Omega = Omega0
J= 3 #dim(Omega)[1]  7
#j= dim(Omega)[2]


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




y_m= vector(mode = "list", length = m)

for(i in 1:m){
for(h in (1:J)){
  
  y_m[[i]][[h]]= rmvnorm(1,rep(0,(dim(K[[h]])[1])),  K[[h]])
  
}
}

#Sum over all the hyperparameters
x_m = lapply(1:m, function(i) Reduce('+', y_m[[i]]))

#Fit a spline
x_m_spl = vector(mode = "list",m)

for(i in 1:m){
  
  x_m_spl[[i]] = splinefun(t, x_m[[i]])
  
}


hat_x_m = vector(mode = "list", m)

for(i in 1:m){
  
  hat_x_m[[i]] = x_m_spl[[i]](t, deriv = 0)   #knots same as obs. points
}


###EMD
dec = vector(mode = "list",m)

for(i in 1:m){
  
   dec[[i]]= emd(hat_x_m[[i]], stoprule = "type2") #, tol = sd(hat_x_m[[i]])*0.1^(20)
   
}



IF= vector(mode = "list", length = m)

for (i in 1: m) {
  for (j in 1: ( dim(dec[[i]]$imf)[2])) {
    IF[[i]][[j]]= ifreq(dec[[i]]$imf[,j], f = 4000)$f[,2]
    
  }
}


IF_hat_x_m<- vector(mode = "list", length = m)

for (i in 1: m) {

  IF_hat_x_m[[i]]= ifreq(hat_x_m[[i]], f = 4000)$f[,2]
    
}




####################
####SPECTROGRAMs####
####################

#m = 1,2,3,4,5
attach(mtcars)
par(mfrow=c(5,3),
    oma = c(6,4,0,0) + 0.1,
    mar = c(0,0,2,1) + 0.1)
s = spectrogram(as.numeric(y_m[[1]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[1]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[1]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')


s = spectrogram(as.numeric(y_m[[2]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[2]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[2]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')


s = spectrogram(as.numeric(y_m[[3]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[3]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[3]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')


s = spectrogram(as.numeric(y_m[[4]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[4]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[4]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')



s = spectrogram(as.numeric(y_m[[5]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[5]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[5]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
title(xlab = "Time(s)",
      ylab = "Frequencies(Hz)",
      outer = TRUE, line = 3)
dev.off()




attach(mtcars)
par(mfrow=c(5,3),
    oma = c(6,4,0,0) + 0.1,
    mar = c(0,0,2,1) + 0.1)
s = spectrogram(as.numeric(y_m[[6]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[6]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[6]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')


s = spectrogram(as.numeric(y_m[[7]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[7]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[7]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')


s = spectrogram(as.numeric(y_m[[8]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[8]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[8]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')


s = spectrogram(as.numeric(y_m[[9]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[9]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[9]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')



s = spectrogram(as.numeric(y_m[[10]][[1]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)

s = spectrogram(as.numeric(y_m[[10]][[2]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(as.numeric(y_m[[10]][[3]]), fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
title(xlab = "Time(s)",
      ylab = "Frequencies(Hz)",
      outer = TRUE, line = 3)
dev.off()







attach(mtcars)
par(mfrow=c(5,4),
    oma = c(6,4,0,0) + 0.1,
    mar = c(0,0,2,1) + 0.1)

s = spectrogram(dec[[1]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)
s = spectrogram(dec[[1]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)
plot(s, yaxt = 'n')
s = spectrogram(dec[[1]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[1]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')


s = spectrogram(dec[[2]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)
s = spectrogram(dec[[2]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)
plot(s, yaxt = 'n')
s = spectrogram(dec[[2]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[2]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')


s = spectrogram(dec[[3]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)
s = spectrogram(dec[[3]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[3]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)  
plot(s, yaxt = 'n')
s = spectrogram(dec[[3]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)  
plot(s, yaxt = 'n')


s = spectrogram(dec[[4]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)  
plot(s)
s = spectrogram(dec[[4]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[4]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[4]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')


s = spectrogram(dec[[5]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)
plot(s)
s = spectrogram(dec[[5]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[5]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[5]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
            timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

title(xlab = "Time(s)",
      ylab = "Frequencies(Hz)",
      outer = TRUE, line = 3)
dev.off()




#m = 6,7,8,9,10
attach(mtcars)
par(mfrow=c(5,4),
    oma = c(6,4,0,0) + 0.1,
    mar = c(0,0,2,1) + 0.1)

s = spectrogram(dec[[6]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)
s = spectrogram(dec[[6]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)
plot(s, yaxt = 'n')
s = spectrogram(dec[[6]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[6]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(dec[[7]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)
s = spectrogram(dec[[7]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)
plot(s, yaxt = 'n')
s = spectrogram(dec[[7]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[7]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(dec[[8]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s)
s = spectrogram(dec[[8]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[8]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)  
plot(s, yaxt = 'n')
s = spectrogram(dec[[8]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)  
plot(s, yaxt = 'n')

s = spectrogram(dec[[9]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)  
plot(s)
s = spectrogram(dec[[9]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[9]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[9]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

s = spectrogram(dec[[10]]$imf[,1], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE)
plot(s)
s = spectrogram(dec[[10]]$imf[,2], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[10]]$imf[,3], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')
s = spectrogram(dec[[10]]$imf[,4], fs =4000 , windowlength = 5, window = 'hamming',
                timestep = -100, dynamicrange = 50, preemphasisf = 50,  maxfreq = 2000, show = FALSE) 
plot(s, yaxt = 'n')

title(xlab = "Time(s)",
      ylab = "Frequencies(Hz)",
      outer = TRUE, line = 3)
dev.off()





################################
### INSTANTANEOUS FREQUENCIES###
################################
#m = 1,2,3,4,5
attach(mtcars)
par(mfrow=c(5,4),
    oma = c(6,4,0,0)+ 0.1, # 
    mar = c(1,1,1,1) + 0.1,
    mgp = c(2,0.5,0)) #
#plot(IF_hat_x_m[[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') 
plot(IF[[1]][[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') #, yaxt = 'n' 
abline(h = median(IF[[1]][[1]]), col = "red")

plot(IF[[1]][[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[1]][[2]]), col = "red")

plot(IF[[1]][[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[1]][[3]]), col = "red")

plot(IF[[1]][[4]], type = 'l', xlab = NA, ylab = NA,xaxt = 'n')
abline(h = median(IF[[1]][[4]]), col = "red")

#plot(IF_hat_x_m[[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') 
plot(IF[[2]][[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') #, yaxt = 'n'
abline(h = median(IF[[2]][[1]]), col = "red")

plot(IF[[2]][[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[2]][[2]]), col = "red")

plot(IF[[2]][[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[2]][[3]]), col = "red")


plot(IF[[2]][[4]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[2]][[4]]), col = "red")

#plot(IF_hat_x_m[[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') 
plot(IF[[3]][[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') #, yaxt = 'n'
abline(h = median(IF[[3]][[1]]), col = "red")

plot(IF[[3]][[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[3]][[2]]), col = "red")

plot(IF[[3]][[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[3]][[3]]), col = "red")

plot(IF[[3]][[4]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[3]][[4]]), col = "red")

#plot(IF_hat_x_m[[4]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') 
plot(IF[[4]][[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') #, yaxt = 'n'
abline(h = median(IF[[4]][[1]]), col = "red")

plot(IF[[4]][[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[4]][[2]]), col = "red")

plot(IF[[4]][[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[4]][[3]]), col = "red")

plot(IF[[4]][[4]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[4]][[4]]), col = "red")

#plot(IF_hat_x_m[[5]], type = 'l', xlab = NA, ylab = NA) 
plot(IF[[5]][[1]], type = 'l', xlab = NA, ylab = NA) #,  yaxt = 'n'
abline(h = median(IF[[5]][[1]]), col = "red")

plot(IF[[5]][[2]], type = 'l', xlab = NA, ylab = NA)
abline(h = median(IF[[5]][[2]]), col = "red")

plot(IF[[5]][[3]], type = 'l', xlab = NA, ylab = NA)
abline(h = median(IF[[5]][[3]]), col = "red")

plot(IF[[5]][[4]], type = 'l', xlab = NA, ylab = NA)
abline(h = median(IF[[5]][[4]]), col = "red")

title(xlab = "Time(s)",
      ylab = "Frequencies(kHz)",
      outer = TRUE, line = 3)
dev.off()



#m = 6,7,8,9,10
attach(mtcars)
par(mfrow=c(5,4),
    oma = c(6,4,0,0)+ 0.1, # 
    mar = c(1,1,1,1) + 0.1,
    mgp = c(2,0.5,0)) #
#plot(IF_hat_x_m[[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') 
plot(IF[[6]][[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') #, yaxt = 'n' 
abline(h = median(IF[[6]][[1]]), col = "red")

plot(IF[[6]][[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[6]][[2]]), col = "red")

plot(IF[[6]][[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[6]][[3]]), col = "red")

plot(IF[[6]][[4]], type = 'l', xlab = NA, ylab = NA,xaxt = 'n')
abline(h = median(IF[[6]][[4]]), col = "red")

#plot(IF_hat_x_m[[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') 
plot(IF[[7]][[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') #, yaxt = 'n'
abline(h = median(IF[[7]][[1]]), col = "red")

plot(IF[[7]][[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[7]][[2]]), col = "red")

plot(IF[[7]][[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[7]][[3]]), col = "red")


plot(IF[[7]][[4]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[7]][[4]]), col = "red")

#plot(IF_hat_x_m[[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') 
plot(IF[[8]][[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') #, yaxt = 'n'
abline(h = median(IF[[8]][[1]]), col = "red")

plot(IF[[8]][[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[8]][[2]]), col = "red")

plot(IF[[8]][[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[8]][[3]]), col = "red")

plot(IF[[8]][[4]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[8]][[4]]), col = "red")

#plot(IF_hat_x_m[[4]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') 
plot(IF[[9]][[1]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n') #, yaxt = 'n'
abline(h = median(IF[[9]][[1]]), col = "red")

plot(IF[[9]][[2]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[9]][[2]]), col = "red")

plot(IF[[9]][[3]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[9]][[3]]), col = "red")

plot(IF[[9]][[4]], type = 'l', xlab = NA, ylab = NA, xaxt = 'n')
abline(h = median(IF[[9]][[4]]), col = "red")

#plot(IF_hat_x_m[[5]], type = 'l', xlab = NA, ylab = NA) 
plot(IF[[10]][[1]], type = 'l', xlab = NA, ylab = NA) #,  yaxt = 'n'
abline(h = median(IF[[10]][[1]]), col = "red")

plot(IF[[10]][[2]], type = 'l', xlab = NA, ylab = NA)
abline(h = median(IF[[10]][[2]]), col = "red")

plot(IF[[10]][[3]], type = 'l', xlab = NA, ylab = NA)
abline(h = median(IF[[10]][[3]]), col = "red")

plot(IF[[10]][[4]], type = 'l', xlab = NA, ylab = NA)
abline(h = median(IF[[10]][[4]]), col = "red")

title(xlab = "Time(s)",
      ylab = "Frequencies(kHz)",
      outer = TRUE, line = 3)
dev.off()



toplot = melt((K[[3]]))
gg = ggplot(toplot, aes(x = Var1, y = Var2,fill = value))
gg = gg + geom_tile()
gg = gg + scale_fill_gradientn(colors = c('lightyellow','blue','red'))
gg = gg + theme_minimal()
gg = gg + theme(panel.border = element_rect(color = 'grey20', fill= NA), panel.grid	= element_blank())
gg = gg + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
gg = gg + theme(legend.position = 'bottom')
gg = gg + labs(x = NULL, y = NULL)
print(gg)

