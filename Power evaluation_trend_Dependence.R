# =============================================================================
#               FONCTION : SIMULATION AND COMPARISON
# =============================================================================


Trend_tests = function(run_num,
                         copula,
                         coefs,
                         Nsim=1000,
                         Nbs=1000,
                         width=10,
                         alpha=0.05
                        )
  {
  
  
  # ETAPE 1 :  PREPARATION
  # ===========================================================================
  
  # Load the necessary packages:
  RequiredPackages <- c("copula","Kendall","VGAM","gtools","openxlsx","resample")

  for (i in RequiredPackages) { 
    if (!require(i, character.only = TRUE)) install.packages(i)
  }
  
  # Set up the progress bar
  wb = txtProgressBar(min=0,max=Nsim+1,style = 3)
  
  
  # Initialization of global variables
  list_of_copulas = list()  # Initialization of a list containing 'couple' objects
  m = length(coefs)         # Number of data to simulate
  T = 1:m                   # Definition of the time vector T
  
  # Initialization of storage for the results of each test
  RATES = rep(0,4)
  
  
  
  # ETAPE 2 : Calculation of copula parameters
  # ===========================================================================
  
  # Definition of which copula to use in data simulation
  switch (copula,
          "clayton" = {
            min_tau = -1
            max_tau = 1
            tau_check_msg = "In this program, the Clayton copula family accepts values of Tau in [-1,1]."
            function_for_copula = function(tau) {claytonCopula(iTau(claytonCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          "gumbel" = {
            min_tau = 0
            max_tau = 1
            tau_check_msg = "In this program, the Gumbel copula family only accepts values of Tau in [0,1)."
            function_for_copula = function(tau) {gumbelCopula(iTau(gumbelCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          "galambos" = {
            min_tau = 0
            max_tau = 1
            tau_check_msg = "In this program, the Galambos copula family only accepts values of Tau in (0,1)."
            function_for_copula = function(tau) {galambosCopula(iTau(galambosCopula(),tau))}
          },
          "joe" = {
            min_tau = 0
            max_tau = 1
            tau_check_msg = "In this program, the Joe copula family only accepts values of Tau in (0,1]."
            function_for_copula = function(tau) {joeCopula(iTau(joeCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          "huslerReiss" = {
            min_tau = 0
            max_tau = 1
            tau_check_msg = "In this program, the Husler-Reiss copula family only accepts values of Tau in (0,1)."
            function_for_copula = function(tau) {huslerReissCopula(iTau(huslerReissCopula(),tau))}
          }, 
          "frank" = {
            min_tau = -1
            max_tau = 1
            tau_check_msg = "In this program, the Frank copula family only accepts values of Tau in (-1,1)."
            function_for_copula = function(tau) {frankCopula(iTau(frankCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
        
  )
  
  # For each element in 'coefs', create a 'copula' object with the corresponding parameter value from 'coefs'
  list_of_copulas = sapply(coefs, function_for_copula)
  
  # ETAPE 3 : DATA SIMULATION
  # ===========================================================================
  # Start of the for-loop
  for (i in 1:Nsim) {
    
    # Update the progress bar
    
    setTxtProgressBar(wb,i)
    
    
    # For each copula parameter, generate 1 observation (x, y)
    rCop = matrix(0,nrow=m,ncol=2)
    for (j in 1:m) {
      rCop[j,] = rCopula(1,list_of_copulas[[j]])
    }
    
    # defini serie X et Y
    X = rCop[,1]
    Y = rCop[,2]

    XY = cbind(X,Y)
    
    
    # ETAPE 4 : CALCULATION OF STATISTICS AND P-VALUES BY METHOD
    # =========================================================================
    
    # Calculation of matrix M, Ck, S, and Cs
    
    # ---------------------
    M = calc_M.mk(XY)
    Ck = calc_C.mk(XY)
    
   
    # ---------------------  
    
    # =========================================================================
   #  METHOD 1: MDT 
    Sn = Sn_method1(X,Y,T,width)
    p.val = bs_pvalue(XY,"MDT",Sn,Nbs,width)
    
    if (p.val < alpha) {
      RATES[1] = RATES[1]+1
    }
    # =========================================================================
    
    # =========================================================================   
    # METHOD 2: MOT
    Sn = Sn_method2(X,Y,T,width)
    p.val = bs_pvalue(XY,"MOT",Sn,Nbs,width)
    
    if (p.val < alpha) {
      RATES[2] = RATES[2]+1
    }
    # =========================================================================
    
 
    
    # =========================================================================    
    # METHOD 4: CIT 
    Sn = Sn_method4(XY,M,Ck)[1]
    
    p.val = MK.CIT_pvalue(Sn,Ck)                  
    if (p.val < alpha) {
      RATES[3] = RATES[3]+1
    }
    # =========================================================================
    
    # =========================================================================    
    # METHOD 5: CET 
    Sn = Sn_method5(XY,M)
    p.val = bs_pvalue(XY,"CET.mk",Sn,Nbs,width)
    
    if (p.val < alpha) {
      RATES[4] = RATES[4]+1
    }
    # =========================================================================
    
    # =========================================================================    
   
  }
  
  PER_RATES = RATES*100/Nsim
  
  
  # ETAPE 5 : SAVE RESULTS
  # ===========================================================================
  
  table = as.data.frame(t(PER_RATES))
  colnames(table) = c("MDT","MOT","CIT","CET ")
  
  if (run_num > 0) {
    line_num = run_num+1

    compare = loadWorkbook("COMPARE.xlsx")
    writeData(compare,sheet="test_power",x=table,colNames = FALSE, rowNames = FALSE,startRow = line_num,startCol = 2)
    
    saveWorkbook(compare,file="COMPARE.xlsx",overwrite = TRUE)
  }
  
  setTxtProgressBar(wb,Nsim+1)
  
  close(wb)
  
  return(table)
}

#coefs=seq(0,1,length.out=m): Data series with trends in dependence.
#coefs=rep(tau,m): Data series without trends (under H0).

