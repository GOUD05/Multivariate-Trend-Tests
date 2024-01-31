# =============================================================================
#               FONCTION : SIMULATION et COMPARAISON
# =============================================================================


compare_tests_dist = function(run_num,
                              copula,
                              coefs,
                              Nsim=1000,
                              Nbs=1000,
                              width=10,
                              alpha=0.05,
                              Ax=0,
                              Bx=0,
                              Ay=0,
                              By=0)
{
  
  
  # ETAPE 1 :  PREPARATION
  # ===========================================================================
  
  # loader des packages necessaires:
  RequiredPackages <- c("copula","Kendall","VGAM","zoo","gtools","openxlsx","resample")
  
  for (i in RequiredPackages) { 
    if (!require(i, character.only = TRUE)) install.packages(i)
  }
  
  # setup pour la bar de progree
  wb = txtProgressBar(min=0,max=Nsim+1,style = 3)
  
  
  # initialization de variables globales
  list_of_copulas = list()  # initalization de liste qui contient les objects "copule"s
  m = length(coefs)         # numero de donnees a simuler 
  T = 1:m                   # definition du vecteur de temps T
  list_of_marginal= list()
  # initialization de storage pour resultats de chaque test
  RATES = rep(0,8)
  
  
  
  # ETAPE 2 : CALCULE DE PARAMETRE DE COPULE
  # ===========================================================================
  
  # definition de quelle copule utilise dans simulation de donnees
  switch (copula,
          "clayton" = {
            #if(any(coefs > 1 || coefs < -1)) {
            #stop("In this program, the Clayton copula family accepts values of Tau in [-1,1].")
            #}
            function_for_copula = function(tau) {claytonCopula(iTau(claytonCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          
          "gumbel" = {
            if(any(coefs >= 1 || coefs < 0)) {
              stop("In this program, the Gumbel copula family only accepts values of Tau in [0,1).")
            }
            function_for_copula = function(tau) {gumbelCopula(iTau(gumbelCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          
          "galambos" = {
            if(any(coefs >= 1 || coefs < 0)) {
              stop("In this program, the Galambos copula family only accepts values of Tau in (0,1).")
            }
            function_for_copula = function(tau) {galambosCopula(iTau(galambosCopula(),tau))}
          },
          
          "joe" = {
            if(any(coefs > 1 || coefs <= 0)) {
              stop("In this program, the Joe copula family only accepts values of Tau in (0,1].")
            }
            function_for_copula = function(tau) {joeCopula(iTau(joeCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          
          "huslerReiss" = {
            if(any(coefs >= 1 || coefs < 0)) {
              stop("In this program, the huslerReissCopula copula family only accepts values of Tau in (0,1).")
            }
            function_for_copula = function(tau) {huslerReissCopula(iTau(huslerReissCopula(),tau))}
          }, 
          
          "frank" = {
            if(any(coefs >= 1 || coefs <= -1)) {
              stop("In this program, the Frank copula family only accepts values of Tau in (-1,1).")
            }
            function_for_copula = function(tau) {frankCopula(iTau(frankCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          
          "amh" = {
            if(any(coefs > 0.3333 || coefs <= -0.1817)) {
              stop("From documentation, the AMH copula family only accepts values of Tau in (-0.1817,0.3333).")
            }
            function_for_copula = function(tau) {amhCopula(iTau(amhCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          
          "tawn" = {
            stop("The Tawn copula family is not yet available.")
            #function_for_copula = function(tau) {tawnCopula(iTau(tawnCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          
          "tev" = {
            stop("The t-EV copula family is not yet available.")
            #function_for_copula = function(tau) {tevCopula(iTau(tevCopula(),tau),dim = 2, use.indepC = "FALSE")}
          },
          
          "normal" = {
            if(any(coefs > 1 || coefs < 1)) {
              stop("In this program, the Normal copula family accepts values of Tau in [-1,1].")
            }
            function_for_copula = function(tau) {normalCopula(iTau(normalCopula(),tau),dim = 2)}
          }
  )
  
  # pour chaque element de "coefs", cree un object "copule" avec le parametre correspondans a la valeur "coefs"
  list_of_copulas = sapply(coefs, function_for_copula)
  
  # ETAPE 3 : SIMULATION de DONNEES
  # ===========================================================================
  # debut de for-loop
  for (i in 1:Nsim) {
    
    # update de la bar de progres 
    setTxtProgressBar(wb,i)
    
    # pour chaque parametre de copule, genere 1 observation (x,y)
    margin_1 = rep(0,m)
    margin_2 = rep(0,m)
    for (j in 1:m) {
      rCop = rCopula(1,list_of_copulas[[j]])
      margin_1[j] = qgev(rCop[1],location = -0.1*j,scale = 1,shape = -0.1)
      margin_2[j] = qgev(rCop[2],location = -0.1*j,scale = 1,shape = -0.1)
    }
    
    # defini serie X et Y
    X = margin_1
    Y = margin_2
    
    # ajouter tendance linear "+ At + B"
    # par default, il n'y a pas de tendance
    #X = X + Ax*T + Bx
    #Y = Y + Ay*T + By
    
    XY = cbind(X,Y)
    
    
    # ETAPE 4 : CALCULE DE STAT ET PVALUE par METHODE
    # =========================================================================
    
    # Calcul de matrice M, Ck, S et Cs
    # METHOD 3/4 on besoin de M, Ck
    # METHOD 5 a besoin seulement de M
    # METHOD 6/7 on besoin de S, Cs
    # METHOD 8 a besoin seulement de S
    
    # WARNING! Calcule de Ck et Cs est tres lent
    
    # ---------------------
    M = calc_M.mk(XY)
    Ck = calc_C.mk(XY)
    
    #S = calc_M.sp(XY)
    #Cs = calc_C.sp(XY)
    # ---------------------  
    
    # =========================================================================
    #  METHOD 1: RWC TEST
    Sn = Sn_method1(X,Y,T,width)
    p.val = bs_pvalue(XY,"RWC",Sn,Nbs,width)
    
    if (p.val < alpha) {
      RATES[1] = RATES[1]+1
    }
    # =========================================================================
    
    # =========================================================================   
    # METHOD 2: RW
    Sn = Sn_method2(X,Y,T,width)
    p.val = bs_pvalue(XY,"rw",Sn,Nbs,width)
    
    if (p.val < alpha) {
      RATES[2] = RATES[2]+1
    }
    # =========================================================================
    
    # =========================================================================   
    # METHOD 3: CST MK
    #Sn = Sn_method3(XY,M,Ck)
    #p.val = bs_pvalue(XY,"CST.mk",Sn,Nbs,width)    # BOOTSTRAP
    #p.val = MK.CST_pvalue(Sn)                      # ASYMP
    
    #if (p.val < alpha) {
    #  RATES[3] = RATES[3]+1
    #}
    # =========================================================================
    
    # =========================================================================    
    # METHOD 4: CIT MK
    Sn = Sn_method4(XY,M,Ck)[1]
    #p.val = bs_pvalue(XY,"CIT.mk",Sn,Nbs,width)   # BOOTSTRAP
    p.val = MK.CIT_pvalue(Sn,Ck)                  # ASYMP
    
    if (p.val < alpha) {
      RATES[4] = RATES[4]+1
    }
    # =========================================================================
    
    # =========================================================================    
    # METHOD 5: CET MK
    #Sn = Sn_method5(XY,M)
    #p.val = bs_pvalue(XY,"CET.mk",Sn,Nbs,width)
    
    #if (p.val < alpha) {
    # RATES[5] = RATES[5]+1
    #}
    # =========================================================================
    
    # =========================================================================    
    # METHOD 6: CST SP
    #Sn = Sn_method6(XY,S,Cs)
    #p.val = bs_pvalue(XY,"CST.sp",Sn,Nbs,width)    # BOOTSTRAP
    #p.val = SP.CST_pvalue(Sn)                      # ASYMP
    
    #if (p.val < alpha) {
    #  RATES[6] = RATES[6]+1
    #}
    # =========================================================================
    
    # =========================================================================    
    # METHOD 7: CIT SP
    #Sn = Sn_method7(XY,S,Cs)
    #p.val = bs_pvalue(XY,"CIT.sp",Sn,Nbs,width)    # BOOTSTRAP
    #p.val = SP.CIT_pvalue(Sn,Cs)                   # ASYMP
    
    #if (p.val < alpha) {
    #  RATES[7] = RATES[7]+1
    #}
    # =========================================================================
    
    # ========================================================================= 
    # METHOD 8: CET SP
    #Sn = Sn_method8(XY,S)
    #p.val = bs_pvalue(XY,"CET.sp",Sn,Nbs,width)
    
    #if (p.val < alpha) {
    #  RATES[8] = RATES[8]+1
    #}
    # =========================================================================
  }
  
  PER_RATES = RATES*100/Nsim
  
  
  # ETAPE 5 : SAUVER RESULTATS
  # ===========================================================================
  
  table = as.data.frame(t(PER_RATES))
  colnames(table) = c("RWC","RW","CST MK","CIT MK","CET MK","CST SP","CIT SP","CET SP")
  
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
