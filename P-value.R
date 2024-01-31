# =============================================================================
#                     FONCTIONS : CALCULATION OF P-VALUES
# =============================================================================


# METHOD 1 : BOOTSTRAP
# =============================================================================

bs_pvalue = function(dat,method,thresh,Nbs,width){
 
  #Length of the series (n)
  n = nrow(dat)
  
# Initialization of TIME series (T)
  T = seq(1:n)
  #Storage of TEST STATISTICS (Sn)
  Sn = rep(0,Nbs)
  
  #Create a matrix (index) of dimensions (n) x (Nbs)
  #Each column of (index) contains the indices for resampling WITH REPLACEMENT.
  index = samp.bootstrap(n,Nbs)
  
  #Using the specified method, run the bootstrap.
  switch(method, 
         "MDT" = {
           for (i in 1:Nbs) {
             dat_sample = dat[index[,i],]
             Sn[i] = Sn_method1(dat_sample[,1],dat_sample[,2],T,width)
           }
         },
         
         "MOT" = {
           for (i in 1:Nbs) {
             dat_sample = dat[index[,i],]
             Sn[i] = Sn_method2(dat_sample[,1],dat_sample[,2],T,width)
           }
         },
         
         
         "CIT.mk" = {
           for (i in 1:Nbs) {
             dat_sample = dat[index[,i],]
             M = calc_M.mk(dat_sample)
             C.mk = calc_C.mk(dat_sample)
             Sn[i] = Sn_method4(dat_sample,M,C.mk)[1]
           }
         },
         
         "CET.mk" = {
           for (i in 1:Nbs) {
             dat_sample = dat[index[,i],]
             M = calc_M.mk(dat_sample)
             Sn[i] = Sn_method5(dat_sample,M)
           }
         }
  )
  return(sum(Sn > rep(thresh,Nbs))/Nbs)
}	


# METHOD 2 : ASYMPTOTIC DISTRIBUTION
# =============================================================================
#This method uses the known distribution of test statistics
#Applies only to the tests MK.CST, MK.CIT

MK.CST_pvalue = function(Sn) {
  pvalue = 2*(1-pnorm(Sn,0,1))
}

MK.CIT_pvalue = function(Sn,C) {
  q = qr(C)$rank
  pvalue = 1 - pchisq(Sn,q)
}

