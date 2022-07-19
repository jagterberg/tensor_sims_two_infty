if(!require(irlba)) {
  install.packages("irlba")
  library(irlba)
}
if(!require(Matrix)) {
  install.packages("Matrix")
  library(Matrix)
}
if(!require(rTensor)) {
  install.packages("rTensor")
  library(rTensor)
}
if(!require(gtools)) {
  install.packages("gtools")
  library(gtools)
}
if(!require(foreach)) {
  install.packages("foreach")
  library(foreach)
}
if(!require(doParallel)) {
  install.packages("doParallel")
  library(doParallel)
}


HOOI_dd <- function(tens,r=c(1,1,1),niter=NULL) {
  to_svd <- k_unfold(tens,1)@data
  to_svd <- to_svd %*% t(to_svd)
  diag(to_svd) <- 0
  init1 <- irlba(to_svd,r[1])$u
  
  to_svd <- k_unfold(tens,2)@data
  to_svd <- to_svd %*% t(to_svd)
  diag(to_svd) <- 0
  init2 <- irlba(to_svd,r[2])$u
  return(HOOI_vanilla(tens=tens,r=r,niter=niter,init1=init1,init2=init2))
}




HOOI_vanilla <- function(tens,r=c(1,1,1),niter=NULL,init1=NULL,init2=NULL) {
  p1 <- dim(tens)[1]
  p2 <- dim(tens)[2]
  p3 <- dim(tens)[3]
  r1 <- r[1]
  r2 <- r[2]
  r3 <- r[3]
  if(is.null(niter)) {
    niter <- max(5,log(max(p1,p2,p3)))
  }
  #get initialization
  if(is.null(init1)) {
    print("obtaining initialization 1 by vanilla (truncated) SVD...")
    init1 <- irlba(k_unfold(tens,1)@data,r1)$u
  } 
  
  if(is.null(init2)) {
    print("obtaining initialization 1 by vanilla (truncated) SVD...")
    init2 <- irlba(k_unfold(tens,2)@data,r2)$u
  }
  
  print("beginning power iteration...")
  
  to_irlba <- ttm( ttm(tens,t(init1),m=1),t(init2),m=2)
  dim3 <- k_unfold(to_irlba,3)@data
  Uhat3 <- svd(dim3,r3)$u
  Uhat1 <- init1
  Uhat2 <- init2
  i <- 1
  while (i < niter) {
    print(paste("iteration",i,"of",niter))
    to_irlba <- ttm(ttm(tens,t(Uhat3),m=3),t(Uhat2),m=2)
    dim1 <- k_unfold(to_irlba,1)@data
    Uhat1 <- svd(dim1)$u[,c(1:r1)]
    
    to_irlba <- ttm(ttm(tens,t(Uhat1),m=1),t(Uhat3),m=3)
    dim2 <- k_unfold(to_irlba,2)@data
    Uhat2 <- svd(dim2)$u[,c(1:r2)]
    
    to_irlba <- ttm(ttm(tens,t(Uhat1),m=1),t(Uhat2),m=2)
    dim3 <- k_unfold(to_irlba,3)@data
    Uhat3 <- svd(dim3)$u[,c(1:r3)]
    i <- i+1
  }
  
  return(list(Uhat1,Uhat2,Uhat3))
}


procrustes <- function(Uhat,U) {
  inner_prod <- t(Uhat) %*% U
  toReturn <- svd(inner_prod)
  
  return(toReturn$u %*% toReturn$v)
}

two_infty <- function(Uhat,U) {
  W <- procrustes(Uhat,U)
  norms <- apply(Uhat %*% W-U,1,function(x){
    sum(x^2)
  })
  return(
    max(sqrt(norms))
  )
}

get_TVZLambda_hat <- function(T_init,uhats) {
  Uhat1 <- uhats[[1]]
  Uhat2 <- uhats[[2]]
  Uhat3 <- uhats[[3]]
  r1 <- dim(Uhat1)[2]
  r2 <- dim(Uhat2)[2]
  r3 <- dim(Uhat3)[2]
  
  #run one iteration:
  to_irlba <- ttm(ttm(T_init,Uhat3%*%t(Uhat3),m=3),Uhat2%*%t(Uhat2),m=2)
  dim1 <- k_unfold(to_irlba,1)@data
  svd1 <- irlba(dim1,r1)
  Uhat1 <- svd1$u
  Vhat1 <- svd1$v
  Lambdahat1 <- svd1$d
  
  #run one iteration:
  to_irlba <- ttm(ttm(T_init,Uhat2%*%t(Uhat2),m=2),Uhat1%*%t(Uhat1),m=1)
  dim2 <- k_unfold(to_irlba,2)@data
  svd2 <- irlba(dim2,r2)
  Uhat2 <- svd2$u
  Vhat2 <- svd2$v
  Lambdahat2 <- svd2$d
  
  #run one iteration:
  to_irlba <- ttm(ttm(T_init,Uhat3%*%t(Uhat3),m=3),Uhat1%*%t(Uhat1),m=1)
  dim3 <- k_unfold(to_irlba,3)@data
  svd3 <- irlba(dim3,r3)
  Uhat3 <- svd3$u
  Vhat3 <- svd3$v
  Lambdahat3 <- svd3$d
  
  V_list <- list(Vhat1,Vhat2,Vhat3)
  Lambdahat_list <- list(Lambdahat1,Lambdahat2,Lambdahat3)
  
  #get estimated T tensor
  T_estimated <- ttm(ttm(ttm(
    T_init,Uhat1%*%t(Uhat1),m=1)
    , Uhat2 %*% t(Uhat2),m=2)
    ,Uhat3 %*% t(Uhat3),m=3)
  
  Z_estimated <- T_init - T_estimated
  
  return(list(T_estimated = T_estimated,
              Z_estimated = Z_estimated,Vs= V_list,
              lambda_hat <- Lambdahat_list))
  
}


get_sijk_hat <- function(i,j,k,T_estimated,Z_estimated,Vs) {
  p1 <- dim(T_estimated)[1]
  p2 <- dim(T_estimated)[2]
  p3 <- dim(T_estimated)[3]
  index_i <- (j-1)*p3 + k
  index_j <- (k-1)*p1 + i
  index_k <- (i-1)*p2 + j
  
  
  Zhat1 <- k_unfold(Z_estimated,1)@data
  Vhat1_proj <- Vs[[1]] %*% Vs[[1]][index_i,]
  Sigma_hat_i <- Zhat1[i,]^2
  sigmahat1 <- sum(
    Vhat1_proj^2 * Sigma_hat_i
  )
  
  Zhat2 <- k_unfold(Z_estimated,2)@data
  Vhat2_proj <- Vs[[2]] %*% Vs[[2]][index_j,]
  Sigma_hat_j <- Zhat2[j,]^2
  sigmahat2 <- sum(
    Vhat2_proj^2 * Sigma_hat_j
  )
  
  Zhat3 <- k_unfold(Z_estimated,3)@data
  Vhat3_proj <- Vs[[3]] %*% Vs[[3]][index_k,]
  Sigma_hat_k <- Zhat2[k,]^2
  sigmahat3 <- sum(
    Vhatk_proj^2 * Sigma_hat_k
  )
  
  shat_ijk <- sqrt(sigmahat1 + sigmahat2 + sigmahat3)
  return(shat_ijk)
  
}

get_Gamma_m_hat <- function(m,mode=1, Z_estimated,Vhats,Lambdas) {
  Z_mode <- k_unfold(Z_estimated,1)@data
  Sigma_hat_m <- Z_mode[m,]^2
  Vhat_mode <- Vhats[[mode]]
  Lambda_mode_inv <- 1/Lambdas[[mode]]
  Gamma_m_hat <- Diagonal(x=Lambda_mode_inv) %*% t(Vhat_mode) %*%
    Diagonal(x=Sigma_hat_m) %*% Vhat_mode %*% Diagonal(x=Lambda_mode_inv)
  
}


