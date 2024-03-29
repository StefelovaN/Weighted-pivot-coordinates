### Weighted pivot coordinates

require("robCompositions")
require("robustHD")

####################################################################################################################################################

weightPivotCoord = function(x, pivotvar = 1, option = "var", method = "classical", 
                            pow = 1, yvar = NULL){
  # INPUT: x - object of class data.frame or matrix; positive values only
  #        pivotvar - pivotal variable; if any other number than 1, the data are resorted in that sense that 
  #                   pivotvar is shifted to the first part
  #        option - option for the choice of weights,
  #                 if option = "var" (default), weights are based on variation matrix elements: (1/t_{1j})^pow,
  #                 if option = "cor", weights are based on correlations between variable specified in yvar and logratios 
  #                              and its distribution:
  #                              |integral_{0}^{r_{j}} f(x) dx|, 
  #                              f(x)... Kernel density estimator for s_{j}; s_{j}=0 if |r_{j}|<cut otherwise s_{j}=r_{j},
  #                              cut = min(#r_{j}=>0/#r_{j}, #r_{j}<0/#r_{j},
  #                              with Gaussian Kernel function and bandwirh h=0.05;
  #        method - method for estimation of variation/corellation,
  #                 if option = "classical" (default), classical estimation is applied,
  #                 if option = "robust", robust estimation is applied;
  #        pow - if option = "var", power pow is applied on unnormalized weights; default is 1;
  #        yvar - if option = "cor", weights are based on correlation between logratios and variable specified in yvar;
  # OUTPUT: WPC - weighted pivot coordinates (matrix with n rows and (D-1) columns);
  #         w - logcontrasts (matrix with D rows and (D-1) columns)
  
  if (option == "var") {
    x1 = cbind(x[,pivotvar], x[,-pivotvar]) 
    if (method == "classical") {
      tj1 = compositions::variation(acomp(x1))[,1]
    } else if (method == "robust") {
      tj1 = robCompositions::variation(acomp(x1))[,1]
    }
    weights2w = 1/tj1[-1]^pow
  } else if (option == "cor") {
    x1 = cbind(x[,pivotvar], x[,-pivotvar]) 
    colnames(x1)[1] = colnames(x)[pivotvar]
    cn = colnames(x1)
    combinations = combn(cn, 2, simplify = FALSE)[1:((ncol(x))-1)]
    codalr = as.data.frame(lapply(combinations, function(j, a) log(a[, j[1]]/a[, j[2]]), a = x1))
    rrj = NULL
    for(j in 1:ncol(codalr)){
      if (method == "classical") {
        rj = cor(yvar, codalr[,j])
      } else if (method == "robust") {
        rj = corHuber(yvar, codalr[,j])
      }
      rrj = c(rrj, rj)
    }
    ssj = rrj
    cut = round(2*min(mean(rrj>=0), mean(rrj<0)), 1)
    ssj[abs(rrj)<cut] = 0
    KernelEstC = function(x){
      x_i = ssj
      h = 0.05
      n = length(x_i)
      return(sum(gauss((x-x_i)/h))/(h*n))
    }
    weights2w = NULL
    D = ncol(x1)
    for(k in 1:(D-1)){
      weights2w = c(weights2w, abs(integrate(f=Vectorize(KernelEstC), lower=0 , upper=rrj[k])$value))
    }
    weights2ww = weights2w  
    for(k in 1:(D-1)){
      vvk = round(weights2w[k], 4)
      vv = round(weights2w[-k], 4)
      if(sum(vv==vvk)>0) weights2ww[k] = weights2ww[k]+abs(rrj[k])/1000
    }
    weights2w = weights2ww
  }
  weights2 = weights2w/sum(weights2w)
  weights1 = c(1, weights2)
  WW = WCoef(weights = weights1)
  WPC = log(x1)%*%WW
  return(list("WPC" = WPC, "w" = WW))
}



############################################################################################################
# the norm of a vector
norm_vec = function(x) sqrt(sum(x^2))

# functions to compute coefficients for construction of weighted coordinates

# w_{1}
vecw1 = function(weights){
  D = length(weights)
  vec = c(1, -weights[-1])
  norm = norm_vec(vec)
  ww1 = vec/norm
  return(ww1)
}

# w_{2}
vecw2 = function(weights){
  D = length(weights)
  vec = c(rep(0, D-3), weights[D-1]-weights[D], weights[D]-weights[D-2], weights[D-2]-weights[D-1])
  norm = norm_vec(vec)
  ww2 = vec/norm
  return(ww2)
}

# w_{p}, p = 3,...,D-2
vecwp = function(weights, p){
  D = length(weights)
  s1 = 0
  for(i in (D-p+1):(D-1)){
    for(j in (i+1):D){
      s1 = s1+weights[i]*weights[j]
    }
  }
  s2 = sum(weights[(D-p+1):D]^2)
  s3 = sum(weights[(D-p+1):D])
  vec = rep(0, D)
  vec[D-p] = 2*s1-(p-1)*s2
  for (k in (D-p+1):D){
    vec[k] = s2-weights[k]^2-(weights[D-p]+weights[k])*(s3-weights[k])+(p-1)*weights[k]*weights[D-p]
  }
  norm = norm_vec(vec)
  wwp = vec/norm
  return(wwp)
}

# w_{D-1}
vecwDm1 = function(weights){
  D = length(weights)
  s1 = 0
  for(i in 2:(D-1)){
    for(j in (i+1):D){
      s1 = s1+weights[i]*weights[j]
    }
  }
  s2 = sum(weights[-1]^2)
  s3 = sum(weights[-1])
  vec = rep(0, D)
  vec[1] = 2*s1-(D-2)*s2
  for (k in 2:D){
    vec[k] = (s2-weights[k]^2)+(1-weights[k])*(s3-weights[k])-(D-2)*weights[k]
  }
  norm = norm_vec(vec)
  wwDm1 = vec/norm
  return(wwDm1)
}

WCoef = function(weights){
  D = length(weights)
  WW = matrix(0, D, D-1)
  WW[, 1] = vecw1(weights = weights)
  WW[, 2] = vecw2(weights = weights)
  for(k in 3:(D-2)){
    WW[, k] = vecwp(weights = weights, p = k)
  }
  WW[, D-1] = vecwDm1(weights = weights)
  return(WW)
}

# Kernel function
gauss = function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
