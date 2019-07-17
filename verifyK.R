  nSeq=c(50,100,150)
  nRep=10^5
  kRWMat=cbind(kRWSeq,kRWSeqN100,kRWSeqN150)
  powerVerifyK = matrix(nrow=length(pSeq),ncol=length(nSeq))
  powerVerifyK2 = powerVerifyK
  
  for(j in 1:length(nSeq))
  {
    marginX = kRWMat[which(pSeq==0.5),j]
    for(i in 1:length(pSeq))
    {
      if(is.finite(kRWSeq[i]))
      {
        tmp = simulateEQRejectionRate(nSeq[j],nSeq[j],pSeq[i],pSeq[i],kRWMat[i,j],alpha,nRep,tableDigits=4)
        powerVerifyK[i,j] = tmp$summaryRslt["EQ_RW_rslt",]
        tmp2 = simulateEQRejectionRate(nSeq[j],nSeq[j],pSeq[i],pSeq[i],marginX,alpha,nRep,tableDigits=4)
        powerVerifyK2[i,j] = tmp2$summaryRslt["EQ_RW_rslt",]
        message("n:",nSeq[j],", p:", pSeq[i], ", power:", powerVerifyK[i,j], ", power2:", powerVerifyK2[i,j])
      }
    }
  }
  save(nSeq,kRWMat,nRep,powerVerifyK,powerVerifyK2,file="powerVerifyK.rda")
