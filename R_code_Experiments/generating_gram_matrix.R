library("EMD")
library("splines")
library("GPfit")
library("seewave")
library("hht")
library(mvtnorm)
library(reshape2)
library(MASS)
library(ggplot2)
#EXPERIMENT 1-----------------------------------------------------------------------------------
#CHECK PERFORMANCES OF KTA WITH KNOWN HYPERPARAMETERS TO TEST IF WE SHOULD EMPLOY SUCH MEASURE AS MOMENT 
#MATCHING TECHNIQUE #REPLACING MLE IN ESTIMATION OF THE HYPERPARAMETER OF EACH GAUSSIAN PROCESS TAKEN INTO ACCOUNT
set.seed(1)

#Define t    (100 time-series) 
n <- 200     #lenght time-series
m <-1000       #number of time-series
t<- vector(mode = "list", m)
t <- lapply(1:m, function(i) t[[i]]<- c(0,sort(runif(n-2,0,10)),10) ) #not equidistant?seq(0,runif(1,0,10), len = n) )    
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
l_list =  c(0.25,1,3,5,0.25,1,3,5)   #2^seq(-2,100)               
nl = length(l_list)

l_alpha_list = expand.grid(l_list[c(1,2,3,4)],c(0.5,3,10))
nla = nrow(l_alpha_list)

# 4) LOCALLY PERIODIC 
l_list =  c(0.25,1,3,5,0.25,1,3,5)   #2^seq(-2,100)               
nl = length(l_list)

l_alpha_list = expand.grid(l_list[c(1,2,3,4)],c(0.5,3,10))
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
#1) KSE kernel
y_KSE<- list()
y_KSE1<- list()

y_KSE<- lapply(1:nl, function(i){
     
   y_KSE1[[i]]<- lapply(1:m,  function(j)rmvnorm(n, mean = rep(0,n), sigma = Sigma_SE[[i]][[j]]))
  
})

rm(y_KSE1)

#2) KPR kernel
y_KPR<- list()
y_KPR1<- list()

y_KPR<- lapply(1:nla, function(i){
  
  y_KPR1[[i]]<- lapply(1:m,  function(j)rmvnorm(n, mean = rep(0,n), sigma = Sigma_PR[[i]][[j]]))
  
})
rm(y_KPR1)

#3) KRQ kernel  
y_KRQ<- list()
y_KRQ1<- list()

y_KRQ<- lapply(1:nla, function(i){
  
  y_KRQ1[[i]]<- lapply(1:m,  function(j)rmvnorm(n, mean = rep(0,n), sigma = Sigma_RQ[[i]][[j]]))
  
})
rm(y_KRQ1)

#4) KlocPR kernel  
y_KlocPR<- list()
y_KlocPR1<- list()

y_KlocPR<- lapply(1:nla, function(i){
  
  y_KlocPR1[[i]]<- lapply(1:m,  function(j)rmvnorm(n, mean = rep(0,n), sigma = Sigma_locPR[[i]][[j]]))
  
})
rm(y_KlocPR1)

#5) KOLS kernel    
y_KOLS<- list()
y_KOLS1<- list()

y_KOLS<- lapply(1:nla, function(i){
  
  y_KOLS1[[i]]<- lapply(1:m,  function(j)rmvnorm(n, mean = rep(0,n), sigma = Sigma_OLS[[i]][[j]]))
  
})
rm(y_KOLS1)

#check that the Frobenius norm is the same for different packages
library(matrixcalc)
frobenius.norm(y_KSE[[1]])
frobenius.norm(t(y_KSE[[1]])* Sigma_SE[[1]])
norm(y_KSE[[1]], type = 'F')
norm( t(y_KSE[[1]]) * Sigma_SE[[1]], type = 'F')


#Computing Kernel Target Alignment KTA -----------------------------------------------------------------------------
#define a function (faster)
#1) KSE kernel - divided according to hyperparameters nl = 6 and #t = m
KTA_KSE<- list()
KTA_KSE1<- list()
KTA_KSE<- lapply(1:nl, function(i)
  KTA_KSE1[[i]]<- lapply(1:m, function (j)
 (norm( t(y_KSE[[i]][[j]]) * Sigma_SE[[i]][[j]], type = 'F'))/(sqrt(norm(y_KSE[[i]][[j]], type = 'F') * norm(Sigma_SE[[i]][[j]], type = 'F')))))
rm(KTA_KSE1)

#2) KPR kernel   
KTA_KPR<- list()
KTA_KPR1<- list()
KTA_KPR<- lapply(1:nla, function(i)
  KTA_KPR1[[i]]<- lapply(1:m, function (j)
    (norm( t(y_KPR[[i]][[j]]) * Sigma_PR[[i]][[j]], type = 'F'))/(sqrt(norm(y_KPR[[i]][[j]], type = 'F') * norm(Sigma_PR[[i]][[j]], type = 'F')))))
rm(KTA_KPR1)

#3) KRQ kernel      
KTA_KRQ<- list()
KTA_KRQ1<- list()
KTA_KRQ<- lapply(1:nla, function(i)
  KTA_KRQ1[[i]]<- lapply(1:m, function (j)
    (norm( t(y_KRQ[[i]][[j]]) * Sigma_RQ[[i]][[j]], type = 'F'))/(sqrt(norm(y_KRQ[[i]][[j]], type = 'F') * norm(Sigma_RQ[[i]][[j]], type = 'F')))))
rm(KTA_KRQ1)

#4) KlocPR kernel   
KTA_KlocPR<- list()
KTA_KlocPR1<- list()
KTA_KlocPR<- lapply(1:nla, function(i)
  KTA_KlocPR1[[i]]<- lapply(1:m, function (j)
    (norm( t(y_KlocPR[[i]][[j]]) * Sigma_locPR[[i]][[j]], type = 'F'))/(sqrt(norm(y_KlocPR[[i]][[j]], type = 'F') * norm(Sigma_locPR[[i]][[j]], type = 'F')))))
rm(KTA_KlocPR1)

#5) KOLS kernel #DOROTA!!    
KTA_KOLS<- list()
KTA_KOLS1<- list()
KTA_KOLS<- lapply(1:nla, function(i)
  KTA_KOLS1[[i]]<- lapply(1:m, function (j)
    (norm( t(y_KOLS[[i]][[j]]) * Sigma_OLS[[i]][[j]], type = 'F'))/(sqrt(norm(y_KOLS[[i]][[j]], type = 'F') * norm(Sigma_OLS[[i]][[j]], type = 'F')))))
rm(KTA_KOLS1)

######---------------------------------------------------------------------------------------
#plots 
#####----------------------------------------------------------------------------------------

#1) KSE Kernel
count_KTA_KSE<- matrix(0, 1, nl)
count<-0

for (i in 1:nl) {
  for(j in 1:m){
    if(KTA_KSE[[i]][[j]] >0.95){
      count<- count +1
      count_KTA_KSE[,i]<-   count} 
    
  }
  
  count<- 0
  
}

df_KSE<- data.frame(l_list, t(count_KTA_KSE))
new_df_KSE<- df_KSE[order(df_KSE$l_list),]

ggplot(data=new_df_KSE, aes(x=factor(l_list), y=t.count_KTA_KSE.)) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Length scale") +
  ylab("#Aligned Kernels > 0.95") +   ggtitle("Square Exponential Kernel") + theme(plot.title = element_text(hjust = 0.5))


#2) KPR Kernel
count_KTA_KPR<- matrix(0, 1, nla)
count<-0

for (i in 1:nla) {
  for(j in 1:m){
    if(KTA_KPR[[i]][[j]] >0.95){
      count<- count +1
      count_KTA_KPR[,i]<-   count} 
    
  }
  
  count<- 0
  
}

df_KPR<- data.frame(l_alpha_list, t(count_KTA_KPR))
new_df_KPR<- df_KPR[order(df_KPR$Var2),]

ggplot(data=new_df_KPR, aes(x=factor(Var1), y=t.count_KTA_KPR., fill = factor(Var2))) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Length scale") +  labs(fill = "Period p") +
  ylab("#Aligned Kernels > 0.95") +   ggtitle("Periodic Kernel") + theme(plot.title = element_text(hjust = 0.5))



#3) KRQ Kernel
count_KTA_KRQ<- matrix(0, 1, nla)
count<-0

for (i in 1:nla) {
  for(j in 1:m){
    if(KTA_KRQ[[i]][[j]] >0.95){
      count<- count +1
      count_KTA_KRQ[,i]<-   count} 
    
  }
  
  count<- 0
  
}

df_KRQ<- data.frame(l_alpha_list, t(count_KTA_KRQ))
new_df_KRQ<- df_KRQ[order(df_KRQ$Var2),]

ggplot(data=new_df_KRQ, aes(x=factor(Var1), y=t.count_KTA_KRQ., fill = factor(Var2))) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Length scale") +  labs(fill = "alfa") +
  ylab("#Aligned Kernels > 0.95") +   ggtitle("Rational Quadratic Kernel") + theme(plot.title = element_text(hjust = 0.5))


#4) KlocPR Kernel
count_KTA_KlocPR<- matrix(0, 1, nla)
count<-0

for (i in 1:nla) {
  for(j in 1:m){
    if(KTA_KlocPR[[i]][[j]] >0.90){
      count<- count +1
      count_KTA_KlocPR[,i]<-   count} 
    
  }
  
  count<- 0
  
}

df_KlocPR<- data.frame(l_alpha_list, t(count_KTA_KlocPR))
new_df_KlocPR<- df_KlocPR[order(df_KlocPR$Var2),]

ggplot(data=new_df_KlocPR, aes(x=factor(Var1), y=t.count_KTA_KlocPR., fill = factor(Var2))) +
  geom_bar(stat="identity", position=position_dodge()) + xlab("Length scale") +  labs(fill = "alfa") +
  ylab("#Aligned Kernels > 0.95") +   ggtitle("Locally Periodic Kernel") + theme(plot.title = element_text(hjust = 0.5))




