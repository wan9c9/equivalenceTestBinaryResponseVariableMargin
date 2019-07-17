# install.packages("xtable") 

library(parallel)
library(xtable)
suppressMessages(library(mvtnorm,quietly = TRUE,warn.conflicts = FALSE))
## Estimate Variance
# always assume muT - muR = eta
# NI: set eta = -abs(eta)
# NS: set eta = abs(eta)
# NI: muT-muR = -eta; NS: muT-muR = eta
RMLE_NI_BR <- function(obsT, obsR, marginX, EPS = 1e-3) 
{
  nT = length(obsT) 
  nR = length(obsR) 
  
  objFun <- function(muR){
    muT = muR+marginX*sqrt(muR*(1-muR))
    
    sum(obsR)*log(muR)+(nR-sum(obsR))*log(1-muR) +
    sum(obsT)*log(muT)+(nT-sum(obsT))*log(1-muT)
  }
  pR <- optimize(objFun, 
                 c(ifelse(marginX>0,EPS,EPS+marginX^2/(1+marginX^2)),
                   ifelse(marginX>0,(1-EPS)/(1+marginX^2),1-EPS)), 
                 maximum = TRUE )
  muR = pR$maximum 
  muT = muR+marginX*sqrt(muR*(1-muR))
  #sigmaR <- sqrt(muR*(1-muR)) 
  #sigmaT <- sqrt(muT*(1-muT)) 
  list(muT=muT,muR=muR)#, sigmaR=sigmaR, sigmaT=sigmaT)
}

var2IIDBP <- function(nT,nR,muT,muR)
{
  muT*(1-muT)/nT + muR*(1-muR)/nR
}

varNIFun <- function(nT,nR,muT,muR,k)
{
  muT*(1-muT)/nT + (sqrt(muR*(1-muR))-k*(0.5-muR))^2/nR
}

varNSFun <- function(nT,nR,muT,muR,k)
{
  muT*(1-muT)/nT + (sqrt(muR*(1-muR))+k*(0.5-muR))^2/nR
}


#' 
#' marginX is the absoluate value
bootstrapPBar <- function(nT,nR,muT,muR,marginX,seed=0,nRep=10^4)
{
  if(!is.null(seed)) set.seed(seed)
  #n = 50; p = 0.5; marginX = 1 ; nRep=10^4
  be = sapply(1:nRep,FUN = function(x){
    xT = rbinom(nT,1,muT); 
    xR = rbinom(nR,1,muR);
    c(unlist(RMLE_NI_BR(xT,xR,-marginX)),
      unlist(RMLE_NI_BR(xT,xR, marginX)) 
    )} 
    )
  val=rowMeans(be) #note that first element is p_R, second is p_T
  list(NI=list(muT=val[1],muR = val[2]),
       NS=list(muT=val[3],muR = val[4])
       )
}

EQ_TRR <- function(nT,nR,muT,muR,marginX,alpha)
{
  delta = (muT-muR)/sqrt(muR*(1-muR))
  pBar = bootstrapPBar(nT,nR,muT,muR,marginX)
  varR = muR*(1-muR)
  varT = muT*(1-muT)
  q = qnorm(1-alpha,lower.tail = TRUE)
  v1 = varNIFun(nT,nR,muT,muR,marginX)
  v1Bar = varNIFun(nT,nR,pBar$NI$muT,pBar$NI$muR,marginX)
  z1 = sqrt(v1Bar/v1)*(q-(delta+marginX)*sqrt(varR/v1Bar))
  
  v2 = varNSFun(nT,nR,muT,muR,marginX)
  v2Bar = varNSFun(nT,nR,pBar$NS$muT,pBar$NS$muR,marginX)
  z2 = sqrt(v2Bar/v2)*(-q-(delta-marginX)*sqrt(varR/v2Bar))
  corr = (varT/nT + (varR - marginX^2*(0.5-muR)^2)/nR)/sqrt(v1*v2)
  corrMat = diag(2)
  corrMat[2,1] = corr
  TRR = pmvnorm(lower=c(z1,-Inf),upper=c(Inf,z2),corr=corrMat)
  rslt=list(delta = delta,
            pBar = pBar,
            z1 = z1,
            z2 = z2,
            corr = corr,
            TRR = TRR
  )
  invisible(rslt)
}


kRWETEq <- function(n,p,alpha,beta,marginX)
{
  rmleEst = bootstrapPBar(n,n,p,p,marginX)
  val=EQ_TRR(n,n,p,p,marginX,alpha)
  val=val$TRR
  #valVec=(c(marginX,rmleEst[2],rmleEst[1],val))
  #names(valVec) =c("marginX","muT","muR","kRWEq"); print(valVec)
  val
}

kRWET <-function(n,p,alpha,beta,bounds=c(0.01,5))
{
  #n = 50; p = 0.9; alpha=0.05;beta=0.1;#marginX = 1 ; nRep=10^4
  tmpFunc <- function(marginX)
  {
    kRWETEq(n,p,alpha,beta,marginX)
  }
  
  #debug
  if(T)
  {
    tmpSeq = seq(bounds[1],bounds[2],by=0.2) 
    valSeq = numeric(length(tmpSeq)) 
    for(i in 1:length(tmpSeq)) 
      valSeq[i] = tmpFunc(tmpSeq[i])
    plot(tmpSeq,valSeq)
  }
  k = uniroot(tmpFunc,bounds)
  k$root
}

if(FALSE)
{
  kRWET(50,0.5,alpha,beta)
  
  pSeq = seq(0.05,0.95,0.05)
  kRWSeq = numeric(length(pSeq)) #kMWO(50,pSeq,alpha,beta)
  for(i in 1:length(pSeq)){
    kRWSeq[i] = tryCatch(kRWET(50,pSeq[i],alpha,beta,c(0.01,ifelse(i==1,2,kRWSeq[i-1]*4))),error=function(e){message("Error in calculating kRW for p",pSeq[i]);NaN})
    message("p:",pSeq[i], ", kRW:",kRWSeq[i])
  } 
  save(pSeq,kMWOSeq,kRWSeq,kRWOSeq,file="calculatedMarginXEQ.rda")
  #EQ_TRR(100,100,0.25,0.5,0.5,0.05)
}

computeEQ_TRR <- function(nT,nR,marginX,alpha)
{
  muRSeq = 0.8
  deltaSeq = seq(-1,1,by=0.25)
  
  TRR = matrix(nrow=length(muRSeq),ncol=length(deltaSeq))
  for(i in 1:length(muRSeq))
  {
    muTSeq = muRSeq[i] + deltaSeq*sqrt(muRSeq[i]*(1-muRSeq[i]))
    for(j in 1:length(muTSeq))
    {
      TRRtmp = EQ_TRR(nT,nR,muTSeq[j],muRSeq[i],marginX,alpha)
      TRR[i,j] = TRRtmp$TRR
      message("muT: ",muTSeq[j], ", muR:", muRSeq[i])
      print(unlist(TRRtmp))
    }
  }
  rownames(TRR)=paste0("muR:",muRSeq)
  colnames(TRR)=paste0("delta:",deltaSeq)
  invisible(TRR)
}

#tst = computeEQ_TRR(50,50,1,0.05)


if(FALSE)
{
  alpha = 0.05
  marginX = 0.26
  muRSeq = seq(0.4,0.8,by=0.2)
  deltaSeq = seq(-1,1,by=0.25)
  
  muTSeq = muRSeq - marginX*sqrt(muRSeq*(1-muRSeq))
  nT = 200
  nR = 200
  
  for(i in 1:length(muRSeq))
  {
    message("muT: ",muTSeq[i], ", muR:", muRSeq[i])
    t = EQ_TRR(nT,nR,muTSeq[i],muRSeq[i],marginX,alpha)
    print(t)
  }
  muRSeq = seq(0.2,0.8,by=0.2)
  muTSeq = muRSeq + marginX*sqrt(muRSeq*(1-muRSeq))
  nT = 200
  nR = 200
  for(i in 1:length(muRSeq))
  {
    message("muT: ",muTSeq[i], ", muR:", muRSeq[i])
    t = EQ_TRR(nT,nR,muTSeq[i],muRSeq[i],marginX,alpha)
    print(t)
  }
  
  TRRMat <- computeEQ_TRR(nT,nR,marginX,alpha)
}


## Wald Test
NI_BR <- function(obsT, obsR, marginX, alpha) 
{
  nT = length(obsT) 
  nR = length(obsR)
  
  smplmuT = mean(obsT)
  smplmuR = mean(obsR) 
  
  absMargin = abs(marginX*sqrt(smplmuR*(1-smplmuR)))
  
  #tstatL_MLE_WO = (smplmuT-smplmuR+absMargin)/sqrt(smplmuT*(1-smplmuT)/nT+smplmuR*(1-smplmuR)/nR) 
  #tstatL_MLE_W = (smplmuT-smplmuR+absMargin)/sqrt(smplmuT*(1-smplmuT)/nT+(sqrt(smplmuR*(1-smplmuR))-marginX*(0.5-smplmuR))^2/nR)
  estL = RMLE_NI_BR(obsT, obsR, -marginX)
  #tstatL_WO = (smplmuT-smplmuR+absMargin)/sqrt(estL$muT*(1-estL$muT)/nT+estL$muR*(1-estL$muR)/nR)
  tstatL_W  = (smplmuT-smplmuR+absMargin)/sqrt(varNIFun(nT,nR,estL$muT,estL$muR,marginX))#(estL$muT*(1-estL$muT)/nT+(sqrt(estL$muR*(1-estL$muR))-marginX*(0.5-estL$muR))^2/nR)
  
  #tstatU_MLE_WO = (smplmuT-smplmuR-absMargin)/sqrt(smplmuT*(1-smplmuT)/nT+smplmuR*(1-smplmuR)/nR) 
  
  #tstatU_MLE_W = (smplmuT-smplmuR-absMargin)/sqrt(smplmuT*(1-smplmuT)/nT+(sqrt(smplmuR*(1-smplmuR))-marginX*(0.5-smplmuR))^2/nR)
  
  estU = RMLE_NI_BR(obsT, obsR, marginX)
  #tstatU_WO = (smplmuT-smplmuR-absMargin)/sqrt(estU$muT*(1-estU$muT)/nT+estU$muR*(1-estU$muR)/nR)
  tstatU_W  = (smplmuT-smplmuR-absMargin)/sqrt(varNSFun(nT,nR,estU$muT,estU$muR,marginX))#sqrt(estU$muT*(1-estU$muT)/nT+(sqrt(estU$muR*(1-estU$muR))-marginX*(0.5+estU$muR))^2/nR)
  
  qntl = qnorm(1-alpha, lower.tail=TRUE)
  
  # non-inferiority test 
  #NI_rslt_MLE_WO = ifelse(!is.finite(tstatL_MLE_WO), 1, ifelse(tstatL_MLE_WO > qntl, 1, 0))
  #NI_rslt_MLE_W = ifelse(tstatL_MLE_W > qntl, 1, 0)
  #NI_rslt_WO = ifelse(tstatL_WO > qntl, 1, 0) 
  NI_rslt_W = ifelse(tstatL_W > qntl, 1, 0)
  
  # non-superiority test
  #NS_rslt_MLE_WO = ifelse(!is.finite(tstatU_MLE_WO), 1, ifelse(tstatU_MLE_WO < -qntl, 1, 0)) 
  #NS_rslt_MLE_W = ifelse(tstatU_MLE_W < -qntl, 1, 0) 
  #NS_rslt_WO = ifelse(tstatU_WO < -qntl, 1, 0) 
  NS_rslt_W = ifelse(tstatU_W < -qntl, 1, 0)
  
  # Equivalence Test
  #EQ_rslt_MLE_WO = ifelse(NI_rslt_MLE_WO & NS_rslt_MLE_WO, 1, 0)
  #EQ_rslt_MLE_W  = ifelse(NI_rslt_MLE_W  & NS_rslt_MLE_W, 1, 0) 
  #EQ_rslt_WO = ifelse(NI_rslt_WO & NS_rslt_WO, 1, 0)
  EQ_rslt_W  = ifelse(NI_rslt_W  & NS_rslt_W, 1, 0)  
  
  c(smplmuT,smplmuR, 
    absMargin, 
    estL$muT,estL$muR,  
    estU$muT,estU$muR, 
    #NI_rslt_MLE_WO, NI_rslt_MLE_W,  NI_rslt_WO, 
    NI_rslt_W, 
    #NS_rslt_MLE_WO, NS_rslt_MLE_W, NS_rslt_WO, 
    NS_rslt_W, 
    #EQ_rslt_MLE_WO, EQ_rslt_MLE_W, EQ_rslt_WO, 
    EQ_rslt_W) 
}

EQ_BR <- function(obsT, obsR, marginX, alpha) 
{
  nT = length(obsT) 
  nR = length(obsR)
  
  smplmuT = mean(obsT)
  smplmuR = mean(obsR) 
  
  absMargin = abs(marginX*sqrt(smplmuR*(1-smplmuR)))
  
  estL = RMLE_NI_BR(obsT, obsR, -marginX)
  estU = RMLE_NI_BR(obsT, obsR, marginX)
  
  qntl = qnorm(1-alpha, lower.tail=TRUE)
  
  tstatL_MWO  = (smplmuT-smplmuR+absMargin)/sqrt(var2IIDBP(nT,nR,smplmuT,smplmuR))
  tstatU_MWO  = (smplmuT-smplmuR-absMargin)/sqrt(var2IIDBP(nT,nR,smplmuT,smplmuR))
  EQ_rslt_MWO  = ifelse((tstatL_MWO > qntl) & (tstatU_MWO < -qntl), 1, 0)  
  
  tstatL_RWO  = (smplmuT-smplmuR+absMargin)/sqrt(var2IIDBP(nT,nR,estL$muT,estL$muR))
  tstatU_RWO  = (smplmuT-smplmuR-absMargin)/sqrt(var2IIDBP(nT,nR,estU$muT,estU$muR))
  EQ_rslt_RWO  = ifelse((tstatL_RWO > qntl) & (tstatU_RWO < -qntl), 1, 0)  
  
  tstatL_RW  = (smplmuT-smplmuR+absMargin)/sqrt(varNIFun(nT,nR,estL$muT,estL$muR,marginX))
  tstatU_RW  = (smplmuT-smplmuR-absMargin)/sqrt(varNSFun(nT,nR,estU$muT,estU$muR,marginX))
  EQ_rslt_RW = ifelse((tstatL_RW > qntl) & (tstatU_RW < -qntl), 1, 0)  
  NI_rslt_RW = ifelse(tstatL_RW > qntl , 1, 0)  
  NS_rslt_RW = ifelse(tstatU_RW < -qntl, 1, 0) 
  c(smplmuT,smplmuR, 
    absMargin, 
    estL$muT,estL$muR,  
    estU$muT,estU$muR, 
    tstatL_RW,tstatU_RW,
    NI_rslt_RW,NS_rslt_RW,
    EQ_rslt_MWO, EQ_rslt_RWO, EQ_rslt_RW
    )
}

simulateEQRejectionRate <- function(nT,nR,muT,muR,marginX,alpha,nRep,tableDigits=4){
  message("Simulating rejection rate for muT=",muT, ", muR=",muR, " with marginX=",marginX, ", nRep=",nRep)
  start_time <- Sys.time()
  
  simEquivTestRep <- function(nT, nR)#, muT, muR, marginX, alpha)
  {
    obsT = rbinom(nT, 1, muT)
    obsR = rbinom(nR, 1, muR)
    testRslt = EQ_BR(obsT, obsR, marginX, alpha)
  }

  testfun <- Vectorize(simEquivTestRep) 
  cl <- makeCluster(detectCores()-1) 
  #get library support needed to run the code
  #clusterEvalQ(cl,library(repsych))
  #put objects in place that might be needed for the code
  clusterExport(cl=cl,varlist=c("RMLE_NI_BR","var2IIDBP","varNIFun","varNSFun", "EQ_BR"))
  clusterExport(cl=cl,varlist=c("nT","nR","muT","muR","marginX","alpha","RMLE_NI_BR","NI_BR","varNIFun","varNSFun","testfun"), envir=environment())
  #... then parallel replicate...
  rslt = parSapply(cl, 1:nRep, function(i,...) { testfun(nT,nR)} ,simplify="array")
  
  #stop the cluster
  stopCluster(cl)
  dimnames(rslt)[[1]] = c( "smplmuT","smplmuR", 
                           "margin", 
                           "muTL","muRL", 
                           "muTU","muRU",
                           "tstatL_RW","tstatU_RW",
                           "NI_RW_rslt","NS_RW_rslt",
                           "EQ_MWO_rslt", 
                           "EQ_RWO_rslt", 
                           "EQ_RW_rslt"
                           )
  dimnames(rslt)[[2]] = paste0("nT",nT,"nR",nR)
  dimnames(rslt)[[3]] = paste0("nRep",1:nRep)
  
  Rslt = apply(rslt, c(1,2), mean)
  summaryStatistics_xtable = xtable(Rslt, 
                      caption = paste0("Summary statistics for equivalence tests for $\\mu_T$=",round(muT,digits = 4)," and $\\mu_R$=",round(muR,digits = 4),"."),  
                      digits=tableDigits,  
                      table.placement ="!h",
                      label=paste0("summary_stat_muT",muT))
  #print(tmp_xtable, scalebox=0.7, caption.placement = "bottom")   
  #print((Rslt))
  average_smplmuT = Rslt["smplmuT",]
  average_smplmuR = Rslt["smplmuR",]
  average_margin = Rslt["margin",]
  muTBar_L = Rslt["muTL",] 
  muRBar_L = Rslt["muRL",] #asymptotic limit of RMLE estimate for muR, below for muT
  
  muTBar_U = Rslt["muTU",]
  muRBar_U = Rslt["muRU",] 
  
  qntl = qnorm(1-alpha, lower.tail=TRUE)
  eqTRRSeq = numeric(length(nT))
  for(i in 1:length(nT))
    eqTRRSeq[i] = (EQ_TRR(nT[i],nR[i],muT,muR,marginX,alpha))$TRR
    
  
  end_time <- Sys.time()
  exec_duration <- end_time - start_time
  invisible(list(indiviudalRslt = rslt,
                 summaryRslt = Rslt,
                 summaryStatistics_xtable = summaryStatistics_xtable, 
                 EQ_RW_TRR = eqTRRSeq, 
                 timeStart = start_time,
                 timeEnd = end_time,
                 timeDuration = exec_duration
                 )
            )
}

if(F)
{
#debug parameters
marginX = 0.262
alpha = 0.05
nRep = 10^5
nT = seq(50,400,by=50)
nR = nT

muR = 0.2
muT = muR - marginX*sqrt(muR*(1-muR))
tableDigits=4

srr = simulateEQRejectionRate(nT,nR,muT,muR,marginX,alpha,nRep)
#print(srr$summaryRslt)
#print(srr$EQ_W_TRR)
print(rbind(srr$summaryRslt,srr$EQ_W_TRR))

}


marginBoundary <- function(p,marginX)
{
  lower <- k^2/(1+k^2)
  upper <- 1/(1+k^2)
}

plotMargin <- function(pSeq,k)
{
  lower = pSeq - k*sqrt(pSeq*(1-pSeq))
  upper = pSeq + k*sqrt(pSeq*(1-pSeq))
  plot(pSeq,lower,type="l",lty=1,ylim=c(-0.5,1.5))
  lines(pSeq,upper,lty=2)
  abline(h=0,col='grey')
  abline(h=1,col="grey")
}

#plotMargin(seq(0.01,0.99,by=0.01),0.262)
