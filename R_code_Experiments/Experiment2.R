library("EMD")
library("splines")
library("GPfit")
library("seewave")
library("hht")
library(mvtnorm)
library(reshape2)
library(MASS)
library(ggplot2)

#EXPERIMENT 2 -----------------------------------------------------------------------------------
#CHECK WITH UNKNOWN PARAMETERS PERFORMANCES OF EMD DETECTING DIFFERENT REGIONS OF THE HYPERPARAMETERS
#OF THE KERNEL FUNCTIONS
set.seed(1)

#Define t    (100 time-series) DOROTA
n <- 200     #lenght time-series
m <-1000       #number of time-series
t<- vector(mode = "list", m)
t <- lapply(1:m, function(i) t[[i]]<- c(0, sort(runif(n-2,0,10)),10))  #seq(0,runif(1,0,10), len = n) )
t_grid<- vector(mode = "list", m)
t_grid = lapply(1:m, function(i) expand.grid(t[[i]], t[[i]]))


#Pick a Kernel:
#KERNEL FUNCTIONS
#---------------------------------------------------------------------------------------------------
# 1) SQUARE EXPONENTIAL
K_SE <- function(x1,x2,l = 1) {exp(- abs(x1-x2)^2 / (2 * l^2) ) }

# 2) PERIODIC WITH SE KERNEL
K_PR <- function(x1,x2,l = 1,p = 1){exp(-2 * sin(pi*abs(x1 - x2) / p)^2 / (l^2) ) }

# 3) RATIONAL QUADRATIC (scale mixture of SE with different l)
K_RQ <- function(x1,x2,l= 1,alpha = 1) { ( 1 + abs(x1-x2)^2/(2 * alpha * l^2)  )^{-alpha} }

# 4) LOCALLY PERIODIC
K_locPR <- function(x1,x2,l = 1,p = 1) {K_PR(x1,x2,l,p) * K_SE(x1,x2,l)}

#5) OLS LINEAR REGRESSION   - DOROTA!
K_OLS <- function(x1,x2,sigma2 = 1){x1 * x2 + sigma2}

#GRID EXPANSIONS ---------------------------------------------------------------------------------
#1) SQUARE EXPONENATIAL
l_list =  c(0.25,1,3,5,7,0.25,1,3,5,7)   #2^seq(-2,100)               
nl = length(l_list)

# 2) PERIODIC WITH SE KERNEL
l_list =  c(0.25,1,3,5,0.25,1,3,5)   #2^seq(-2,100)               
nl = length(l_list)

l_alpha_list = expand.grid(l_list[c(1,2,3,4)],c(0.5,3,10))
nla = nrow(l_alpha_list)

# 3) RATIONAL QUADRATIC 
l_list =  c(0.25,1,3,0.25,1,3)   #2^seq(-2,100)               
nl = length(l_list)

l_alpha_list = expand.grid(l_list[c(1,2,3)],c(0.5,7))
nla = nrow(l_alpha_list)

# 4) LOCALLY PERIODIC 
l_list =  c(0.25,1,3,0.25,1,3)   #2^seq(-2,100)               
nl = length(l_list)

l_alpha_list = expand.grid(l_list[c(1,2,3)],c(0.5,7))
nla = nrow(l_alpha_list)

# 5) OLS LINEAR REGRESSION
sigma2_list =  c(0.25,1,3,0.25,1,3)   #2^seq(-2,100)               
nsigma2 = length(sigma2_list)

#Kernel evaluation------------------------------------------------------------------------------
###1) SQUARE EXPONENTIAL

function_Sigma_SE<- function(l,j)lapply(1:m, function(j) (matrix( apply(expand.grid(t[[j]],t[[j]]),1, 
                                                                        function(s) K_SE(x1=s[1],x2=s[2],l = l)),nrow = n)))
#list divided per hyperparameters, i.e. nl = 6 times m (# of time series)
Sigma_SE<- lapply(l_list, function(l,j)function_Sigma_SE(l,j))

###2) PERIODIC WITH SE KERNEL
function_Sigma_PR = function(l,a,j)lapply(1:m, function(j) matrix( apply(expand.grid(t[[j]],t[[j]]),1, 
                                                                         function(s) K_PR(x1=s[1],x2=s[2],l = l, p = a)),nrow = n))

Sigma_PR = lapply(1:nla, function(i) function_Sigma_PR(l_alpha_list[i,1],l_alpha_list[i,2] /2,j) )

### 3) RATIONAL QUADRATIC 
function_Sigma_RQ = function(l,a,j)lapply(1:m, function(j) matrix( apply(expand.grid(t[[j]],t[[j]]),1, 
                                                                         function(s) K_RQ(x1=s[1],x2=s[2],l = l, alpha = a)),nrow = n) )


Sigma_RQ = lapply(1:nla, function(i) function_Sigma_RQ(l_alpha_list[i,1],l_alpha_list[i,2],j ) )

### 4) LOCALLY PERIODIC
function_Sigma_locPR = function(l,a,j) lapply(1:m, function(j) matrix( apply(expand.grid(t[[j]],t[[j]]),1, 
                                                                             function(s) K_locPR(x1=s[1],x2=s[2],l = l, p = a)),nrow = n))


Sigma_locPR = lapply(1:nla, function(i) function_Sigma_locPR(l_alpha_list[i,1],l_alpha_list[i,2]/2 ,j ) )

### 5) OLS LINEAR REGRESSION
function_Sigma_OLS = function(sigma,j) lapply(1:m, function(j) matrix( apply(expand.grid(t[[j]],t[[j]]),1,
                                                                             function(s) K_OLS(x1=s[1],x2=s[2],sigma2 = sigma)),nrow = n))

Sigma_OLS = lapply(sigma2_list, function(sigma) function_Sigma_OLS (sigma,j) )


#Generate multivariate normal distribution ------------------------------------------------------------------------
#Assumption: each IMF is distributed according to a GP(0,K(t,t))
#################################################
#1) KSE kernel
y_KSE<- list()
y_KSE1<- list()

y_KSE<- lapply(1:nl, function(i){
  
  y_KSE1[[i]]<- lapply(1:m,  function(j)rmvnorm(1, mean = rep(0,n), sigma = Sigma_SE[[i]][[j]])) 
  
})

rm(y_KSE1)

#fit a spline
#new_t<- lapply(1:m, function(i) sort(sample(t[[i]], size = n/2))) # interpolation not knots selection
#index<- lapply(1:m, function(i) match(new_t[[i]], t[[i]]))             

#y_KSE1<- list()
#y_KSE_new<- lapply(1:nl, function(i){
  
#  y_KSE1[[i]]<- lapply(1:m, function(j) y_KSE[[i]][[j]][index[[j]]]   )
  
#})
#rm(y_KSE1)

y_KSE_spl<- vector(mode="list", nl)
y_KSE_spl1<- list()
y_KSE_spl<- lapply(1:nl, function(i){
  
    y_KSE_spl1[[i]]<- lapply(1:m, function(j) splinefun(t[[j]], y_KSE[[i]][[j]]))  #splinefun(new_t[[j]], y_KSE_new[[i]][[j]]))
  
})
rm(y_KSE_spl1)

y_KSE_st1<-list()
y_KSE_st<- vector(mode = "list", nl)

y_KSE_st<- lapply(1:nl, function(i){
  
      y_KSE_st1[[i]]<- lapply(1:m, function(j)y_KSE_spl[[i]][[j]](t[[j]], deriv = 0) )  #(new_t[[j]], deriv = 0) )   
  
})

rm(y_KSE_st1)

##########################################
#2) KPR kernel
y_KPR<- list()
y_KPR1<- list()

y_KPR<- lapply(1:nl, function(i){
  
  y_KPR1[[i]]<- lapply(1:m,  function(j)rmvnorm(1, mean = rep(0,n), sigma = Sigma_PR[[i]][[j]])) 
  
})

rm(y_KPR1)

#fit a spline
#new_t<- lapply(1:m, function(i) sort(sample(t[[i]], size = n/2))) # interpolation not knots selection
#index<- lapply(1:m, function(i) match(new_t[[i]], t[[i]]))             

#y_KPR1<- list()
#y_KPR_new<- lapply(1:nl, function(i){

#  y_KPR1[[i]]<- lapply(1:m, function(j) y_KPR[[i]][[j]][index[[j]]]   )

#})
#rm(y_KPR1)

y_KPR_spl<- vector(mode="list", nl)
y_KPR_spl1<- list()
y_KPR_spl<- lapply(1:nl, function(i){
  
  y_KPR_spl1[[i]]<- lapply(1:m, function(j) splinefun(t[[j]], y_KPR[[i]][[j]]))  #splinefun(new_t[[j]], y_KPR_new[[i]][[j]]))
  
})

rm(y_KPR_spl1)

y_KPR_st1<-list()
y_KPR_st<- vector(mode = "list", nl)

y_KPR_st<- lapply(1:nl, function(i){
  
  y_KPR_st1[[i]]<- lapply(1:m, function(j)y_KPR_spl[[i]][[j]](t[[j]], deriv = 0) )  #(new_t[[j]], deriv = 0) )   
  
})

rm(y_KPR_st1)


#ONE SYNTHETIC OBSERVATION SET
#1) KSE y

y_final_KSE11<- list()
y_final_KSE1<- vector(mode = "list", nl)

y_final_KSE1<- lapply(1:m, function(i){
  
    y_final_KSE11[[i]]<- sapply(1:nl, function(j) list(y_KSE_st[[j]][[i]]))
  
})

rm(y_final_KSE11)

y_final_KSE<- lapply(1:m, function(i) Reduce('+', y_final_KSE1[[i]]))
rm(y_final_KSE1)

#2) KPR y

y_final_KPR11<- list()
y_final_KPR1<- vector(mode = "list", nl)

y_final_KPR1<- lapply(1:m, function(i){
  
  y_final_KPR11[[i]]<- sapply(1:nl, function(j) list(y_KPR_st[[j]][[i]]))
  
})

rm(y_final_KPR11)

y_final_KPR<- lapply(1:m, function(i) Reduce('+', y_final_KPR1[[i]]))
rm(y_final_KPR1)

#APPLY THE EMD
#1) KSE DECOMPOSITION BASIS
dec_KSE<- list()
dec_KSE<- lapply(1:m, function(i)emd(y_final_KSE[[i]], t[[i]]) )  #new_t


IMF_KSE<-list()

for(i in(1:m)){
  
  IMF_KSE[[i]]<- dec_KSE[[i]]$imf[,1:dim(dec_KSE[[i]]$imf)[2]]
  
  
}

#1) KPR DECOMPOSITION BASIS
dec_KPR<- list()
dec_KPR<- lapply(1:m, function(i)emd(y_final_KPR[[i]], t[[i]]) )  #new_t


IMF_KPR<-list()

for(i in(1:m)){
  
  IMF_KPR[[i]]<- dec_KPR[[i]]$imf[,1:dim(dec_KPR[[i]]$imf)[2]]
  
  
}
#COMPUTE THE INSTANTANEOUS FREQUENCY  (check which package is better - ifreq(seewave) or InstantaneousF(hht))
#1)KSE IF
IF_KSE<- list()
IF_KSE1<- list()

IF_KSE<- lapply(1:m, function(i){
  
  IF_KSE1[[i]]<- sapply(1:dim(IMF_KSE[[i]])[2], function(j) ifreq( IMF_KSE[[i]][,j],f = 4)$f[,2]) #only IF 
                                                            #how to set f = 4 since the minum period is 0.25s
                                                            #corresponding to 
})

rm(IF_KSE1)


#2)KPR IF
IF_KPR<- list()
IF_KPR1<- list()

IF_KPR<- lapply(1:m, function(i){
  
  IF_KPR1[[i]]<- sapply(1:dim(IMF_KPR[[i]])[2], function(j) ifreq( IMF_KPR[[i]][,j],f = 4)$f[,2]) #only IF 
})

rm(IF_KPR1)



#COMPUTE RANGES FOR INSTANTANEOUS FREQUENCIES to relate it to the periods

range_IF<- list()
range_IF1<- list()


range_IF<- lapply(1:m, function(i){
  
   range_IF1[[i]]<- lapply(1:dim(IF_KPR[[i]])[2], function(j) range(IF_KPR[[i]][,j]))
  
})

rm(range_IF1)


upper_rangeIF<- list()
upper_rangeIF1<- list()

upper_rangeIF<- lapply(1:m, function(i){
  
   upper_rangeIF1<- lapply(1:length(range_IF[[i]]), function(j) max(range_IF[[i]][[j]] *1000)) #Hz not kHz
  
})

rm(upper_rangeIF1)


lower_rangeIF<- list()
lower_rangeIF1<- list()

lower_rangeIF<- lapply(1:m, function(i){
  
  lower_rangeIF1<- lapply(1:length(range_IF[[i]]), function(j) min(range_IF[[i]][[j]] *1000))  #Hz not kHz
  
})

rm(lower_rangeIF1)


count_IF1<- 0

for (i in 1:m) {
  
  if(upper_rangeIF[[i]][1] <  4 ){ count_IF1 = count_IF1+1}
  
}
 
count_IF2<- 0

for (i in 1:m) {
  
  if(upper_rangeIF[[i]][2] <  4  & upper_rangeIF[[i]][2] >  0.6){ count_IF2 = count_IF2+1}
  
}



##############PLOT ----------------------------------------
#1)KSE ----------------------------------------------------
par(mfrow = c(1, 1))
plot(NA, ylim = c(0, 2), xlim = c(0, 200)) #better view over the lines

cl <- rainbow(m)
for(i in 1:m){
  lines(IF_KSE[[i]][,1]*1000,  col = cl[i], lwd=1.5) 
}


par(mfrow = c(1, 1))
plot(NA, ylim = c(0, 2), xlim = c(0, 200)) #better view over the lines

cl <- rainbow(m)
for(i in 1:m){
  lines(IF_KSE[[i]][,2]*1000,  col = cl[i], lwd=1.5) 
}


#1)KPR ----------------------------------------------------
par(mfrow = c(1, 1))
plot(NA, ylim = c(0, 2), xlim = c(0, 200),  ylab = "IF1", xlab = "time") 

cl <- rainbow(m)
for(i in 1:m){
  lines(IF_KPR[[i]][,1]*1000,  col = cl[i], lwd=1.5) 
}


par(mfrow = c(1, 1))
plot(NA, ylim = c(0, 2), xlim = c(0, 200), ylab = "IF2", xlab = "time") 

cl <- rainbow(m)
for(i in 1:m){
  lines(IF_KPR[[i]][,2]*1000,  col = cl[i], lwd=1.5) 
}


par(mfrow = c(1, 1))
plot(NA, ylim = c(0, 2), xlim = c(0, 200),  ylab = "IF3", xlab = "time") 

cl <- rainbow(m)
for(i in 1:m){
  lines(IF_KPR[[i]][,3]*1000,  col = cl[i], lwd=1.5) 
}
abline( h = 0.667, col = "black", lwd = 0.5)

par(mfrow = c(1, 1))
plot(NA, ylim = c(0, 2), xlim = c(0, 200), ylab = "IF4", xlab = "time") 

cl <- rainbow(m)
for(i in 1:m){
  lines(IF_KPR[[i]][,4]*1000,  col = cl[i], lwd=1.5) 
}
abline( h = 0.04, col = "black", lwd = 0.5)


par(mfrow = c(1, 1))
plot(NA, ylim = c(0, 2), xlim = c(0, 200), ylab = "IMFs", xlab = "time") 

cl <- rainbow(m)
for(i in 1:m){
  lines(IF_KPR[[i]][,5]*1000,  col = cl[i], lwd=1.5) 
}

# #COMPUTE THE FFT
# #1) KSE FFT
# FFT_KSE<- lapply(1:nl, function(i) fft( dec_KSE$imf[,i]))
# 
# #2) KPR FFT
# FFT_KPR<- lapply(1:dim(dec_KPR$imf)[2], function(i) fft( dec_KPR$imf[,i]))








  
  