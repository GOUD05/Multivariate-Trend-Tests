# =============================================================================
#                 FONCTIONS : CALCULE DE STATISTIQUE DE TESTE
# =============================================================================


# METHOD 1 : MDT TEST
# =============================================================================

Sn_method1= function(X,Y,T,width) {
  
  T2 = 1:(length(X)-width+1)

  mat = cbind(X,Y)
  roll = rollapply(mat, width = width,FUN = function(x){kendall.tau(x[, 1], x[, 2])}, 
                   by.column = FALSE, align = "right")
  dep1 = kendall.tau(as.numeric(roll),T2)
  return(abs(dep1))
}	


# METHOD 2 : MOT
# =============================================================================

Sn_method2 = function(X,Y,T,width) {
  
  # declaration de variable
  T2 = 1:(length(X)-width+1)
  
  # obtien coefficient de correlation entre 'x' et 'y' sur different sous-ensembles de largeur 'width'
  # XY est une serie de longeur: length(X) - width + 1
  mat = cbind(X,Y)

  roll = rollapply(mat, width = width,FUN = function(x){kendall.tau(x[, 1], x[, 2])}, 
                 by.column = FALSE, align = "right")
  

  # obtien trois mesure de Kendall' tau
  dep1 = (kendall.tau(as.numeric(roll),T2))^2
  dep2 = (kendall.tau(X,T))^2
  dep3 = (kendall.tau(Y,T))^2
  
  # retour la somme div. par 3
  return(abs(1/3*(sum(dep1,dep2,dep3))))
}	



# METHOD 4 : CIT 
# =============================================================================
Sn_method4 = function(data,M,C) {
  D = M %*% solve(C) %*% t(M)
  return(D)
}



# METHOD 5 : CET 
# =============================================================================
Sn_method5 = function(data,M) {
  L = M %*% t(M) 
  return(L)
}


# CALCULE DE MATRICE POUR MK
# =============================================================================

calc_C.mk = function(data) {
  
  size = dim(data)
  
  # Calcul de t(uv)
  t = matrix(nrow = size[2], ncol=size[2], byrow=TRUE)
  
  for (v in 1:size[2]) {
    for (u in 1:size[2]) {
      t[u,v] = 0
      for (j in 2:size[1]) {
        for (i in 1:(j-1)) {
          t[u,v] = t[u,v] + sign((data[j,u]-data[i,u])*(data[j,v]-data[i,v]))
        }
      }
    }
  }
  
  
  # Calcul de r(uv)
  r = matrix(nrow = size[2], ncol=size[2], byrow=TRUE)
  
  for (v in 1:size[2]) {
    for (u in 1:size[2]) {
      r[u,v] = 0
      for ( h in 1:size[1]) {
        for (j in 1:size[1]) {
          for (i in 1:size[1]) {
            r[u,v] = r[u,v] + sign((data[h,u] - data[j,u])*(data[h,v]-data[i,v]))
          }
        }
      }
    }
  }
  
  # Calcul de C(uv)
  C = t/3 + r/3
  for (i in 1:size[2]) {
    C[i,i] = (size[1]*(size[1]-1)*(2*size[1]+5))/18
  }
  
  return(C)
}

calc_M.mk = function(data) {
  
  size = dim(data)
  
  M = matrix(nrow=1, ncol=size[2], byrow=TRUE)
  
  for (u in 1:size[2]) {
    M[u] = 0
    for (j in 2:size[1]) {
      for (i in 1:(j-1)) {
        M[u] = M[u] + sign(data[j,u]-data[i,u])
      }
    }
  }
  
  return(M) 
}
