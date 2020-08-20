# File:         learn_diagnostics 
# Description:  MCMC 
#               for Project 3, Artificial Intelligence 2019, UU
# Author:       Rex Ruan

# Install the package
# install.packages("Diagnostics_1.2.0.tar", repos = NULL, type="source")

library(Diagnostics)

learn <- function(data) {
  # Pn, Te, VTB, TB, Sm, LC, Br, XR, Dy
  
  #Pn, VTB, Sm
  pn = catePro(data$Pn)
  vtb = catePro(data$VTB)
  sm = catePro(data$Sm)
  #TB, LC, Br, Te
  tb = conditwoPro(data$TB,data$VTB)
  lc = conditwoPro(data$LC,data$Sm)
  br = conditwoPro(data$Br,data$Sm)
  te = conditwoPro(data$Te,data$Pn,continous = T)
  #Dy
  dy = condithreePro(data$Dy,data$LC,data$Br)
  #XR
  xr = condifourPro(data$XR,data$Pn, data$TB, data$LC)
  
  return(list(pn=pn,vtb=vtb,sm=sm,tb=tb,lc=lc,br=br,te=te,dy=dy,xr=xr))
}

diagnose <- function(net,cases) {
  randomNumers <- matrix(runif(400000),nrow=10,ncol=40000)
  estimates <- matrix(nrow=10,ncol=4)
  for(i in 1:10){
    estimates[i,] = getProbs(net,cases[i,],randomNumers[i,])
  }
  return(estimates)
}

catePro <- function(cate){
  m <- matrix(ncol=2)
  m[1,2] <- sum(cate)/length(cate)
  m[1,1] <- 1-m[1,2]
  return(m)
}

conditwoPro <- function(B,A,continous=F) {
  m <- matrix(ncol=3, nrow=2)
  m[,1] <- c(0,1)
  if(continous==F){
  m[1,3] <- sum(B[which(A==0)])/length(which(A==0))
  m[1,2] <- 1 - m[1,3]
  m[2,3] <- sum(B[which(A==1)])/length(which(A==1))
  m[2,2] <- 1 - m[2,3]} else {
    tem <- B[which(A==0)]
    m[1,2:3] <- c(mean(tem),sd(tem))
    tem <- B[which(A==1)]
    m[2,2:3] <- c(mean(tem),sd(tem))
  }
  return(m)
}

condithreePro <- function(D,B,C) {
  m <- matrix(ncol=4,nrow=4)
  m[,1] <- c(0,0,1,1)
  m[,2] <- c(0,1,0,1)
  for(i in 1:4){
    index <- intersect(which(B==m[i,1]),which(C==m[i,2]))
    m[i,4] <- sum(D[index])/length(index)
    m[i,3] <- 1 - m[i,4]
  }
  return(m)
}

condifourPro <- function(D,A,B,C) {
  m <- matrix(ncol=5,nrow=8)
  m[,1] <- c(0,0,0,0,1,1,1,1)
  m[,2] <- c(0,0,1,1,0,0,1,1)
  m[,3] <- c(0,1,0,1,0,1,0,1)
  for(i in 1:8){
    index <- intersect(intersect(which(A==m[i,1]),which(B==m[i,2])),which(C==m[i,3]))
    m[i,5] <- sum(D[index])/length(index)
    m[i,4] <- 1 - m[i,5]
  }
  return(m)
}

candFuntion <- function(l,x) {
  if(x=='pn'){
    l$rand <- 1
    if(l$pn==0){l$pn=1}else{l$pn=0}}
  else if(x=='tb'){
    l$rand <- 2
    if(l$tb==0){l$tb=1}else{l$tb=0}}
  else if(x=='lc'){
    l$rand <- 3
    if(l$lc==0){l$lc=1}else{l$lc=0}}
  else if(x=='br'){
    l$rand <- 4
    if(l$br==0){l$br=1}else{l$br=0}}
  return(l)
}

computeProd <- function(net,values) {
  return(prod(
    net$pn[values$pn+1],dnorm(values$te,net$te[values$pn+1,2],net$te[values$pn+1,3]),
    net$vtb[values$vtb+1], net$tb[which(net$tb[,1]==values$vtb),values$tb+2],
    net$sm[values$sm+1], net$lc[which(net$lc[,1]==values$sm),values$lc+2],
    net$br[which(net$br[,1]==values$sm),values$br+2],
    net$xr[intersect(intersect(which(net$xr[,1]==values$pn),which(net$xr[,2]==values$tb)),which(net$xr[,3]==values$lc)),values$xr+4],
    net$dy[intersect(which(net$dy[,1]==values$lc),which(net$dy[,2]==values$br)),values$dy+3]
  ))
}

assProp <- function(net,values,x,rands){
  pOld = computeProd(net,values)
  newValues = candFuntion(values,x)
  pNew = computeProd(net,newValues)
  if(pNew < pOld){
    if(pNew/pOld>rands[newValues$rand]){values<-newValues}
  } else {
    values <- newValues
  }

  return(values)
}  

#order:             pn, te, vtb,  tb, sm, lc, br, xr, dy
#case:                  te, vtb,      sm,         xr, dy
#conditional probs: pn,           tb,     lc, br
#pn, te|pn, vtb, tb|vtb, sm, lc|sm, br|sm, xr|(pn,tb,lc), dy|(lc,br)
generateSample <- function(net,case,rands,initial=F,values=NA){
  if(initial==T){
    assignedValues = sample(0:1,4,replace=T) 
    values <- list(pn=assignedValues[1],te=case$Te,vtb=case$VTB,
                 tb=assignedValues[2],sm=case$Sm,lc=assignedValues[3],
                 br=assignedValues[4],xr=case$XR,dy=case$Dy)}
  for(x in c('pn','tb','lc','br')){
    values = assProp(net,values,x,rands)
  }
  return(values)
}

getProbs <- function(net,case,randomNumbers) {
  randBurnMatrix <- matrix(randomNumbers[1:4000],nrow=1000,ncol=4)
  randSampleMatrix <- matrix(randomNumbers[4001:40000],nrow=9000,ncol=4)
  m <- matrix(ncol=9,nrow=9000)
  #initialization
  values = generateSample(net,case,randBurnMatrix[1,],initial = T)
  #burning period
  for(i in 2:1000){
    values = generateSample(net,list(te=case$Te, vtb=case$VTB,
                                     sm=case$Sm, xr=case$XR,
                                     dy=case$Dy),
                            randBurnMatrix[i,],values=values)
  }
  #sampling period
  for(i in 1:9000){
    values = generateSample(net,list(te=case$Te, vtb=case$VTB,
                                     sm=case$Sm, xr=case$XR,
                                     dy=case$Dy),
                            randSampleMatrix[i,],values=values)
    
    m[i,]<-c(values$pn,values$te,values$vtb,values$tb,values$sm,
             values$lc,values$br,values$xr,values$dy)
  }
  pn = sum(m[,1])/9000
  tb = sum(m[,4])/9000
  lc = sum(m[,6])/9000
  br = sum(m[,7])/9000
  return(c(pn,tb,lc,br))
}


