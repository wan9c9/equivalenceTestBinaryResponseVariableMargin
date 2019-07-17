alpha=0.05;beta=0.10
  pSeq = seq(0.05,0.95,0.05)
  
  kRWSeq = numeric(length(pSeq)) 
  for(i in 1:length(pSeq)){
    kRWSeq[i] = tryCatch(kRWET(50,pSeq[i],alpha,beta,c(0.01,ifelse(i==1,4,kRWSeq[i-1]*4))),error=function(e){message("Error in calculating kRW for p",pSeq[i]);NaN})
    message("p:",pSeq[i], ", kRW:",kRWSeq[i])
  } 
  
  kRWSeqN100 = numeric(length(pSeq)) 
  for(i in 1:length(pSeq)){
    kRWSeqN100[i] = tryCatch(kRWET(100,pSeq[i],alpha,beta,c(0.01,ifelse(i==1,2,kRWSeqN100[i-1]*4))),error=function(e){message("Error in calculating kRW for p",pSeq[i]);NaN})
    message("p:",pSeq[i], ", kRW:",kRWSeqN100[i])
  } 
  
  kRWSeqN150 = numeric(length(pSeq))
  for(i in 1:length(pSeq)){
    kRWSeqN150[i] = tryCatch(kRWET(150,pSeq[i],alpha,beta,c(0.01,ifelse(i==1,1,kRWSeqN150[i-1]*2))),error=function(e){message("Error in calculating kRW for p",pSeq[i]);NaN})
    message("p:",pSeq[i], ", kRW:",kRWSeqN150[i])
  } 
  
  save(alpha,beta,pSeq,kRWSeq,kRWSeqN100,kRWSeqN150,file="calculatedMarginXEQ.rda")
