set.seed(5152022)
numcores<- detectCores()
cl <- makeCluster(4,type = "FORK")
registerDoParallel(cl,numcores)

two_infty_sim <- function(p,r,C,sigmas,ntrials) {
  
  #generate pure nodes
  purenodes <- c(1,rep(0,r-1))
  for(i in c(2:r)) {
    newnode <- c( rep(0,(i-1)),1,rep(0,(r-i)))
    purenodes <- rbind(purenodes,newnode)
  }
  #generate memberships
  Pi1 <- rbind(purenodes,rdirichlet(p-r,rep(1,r)))
  Pi2 <- rbind(purenodes,rdirichlet(p-r,rep(1,r)))
  Pi3 <- rbind(purenodes,rdirichlet(p-r,rep(1,r)))
  
  #true tensor: underlying r x r x r rank tensor
  T_true <- ttm(ttm(ttm(C,Pi1,m=1),Pi2,m=2),Pi3,m=3)
  T1 <- k_unfold(T_true,1)@data
  Utrue1 <- irlba(T1,3,3)$u
  print("got here 2!")
  # for each noise level
  toreturn <- foreach(sigmavals= c(1:length(sigmas)),.combine=rbind) %dopar% {
    sds <- runif(p^3,0,sigmas[sigmavals])
    vals <- foreach(trial = c(1:ntrials),.combine=cbind) %dopar% {
      Z <- as.tensor(array(rnorm(p^3,sd=sds),dim=c(p,p,p)))
      T_obs <- T_true + Z
      Uhats <- HOOI_dd(T_obs,r=c(3,3,3))
      Uhat1 <- Uhats[[1]]
      val <- two_infty(Uhat1,Utrue1)
      return(val)
    }
    return(vals)
  }
  rownames(toreturn) <- sigmas
  colnames(toreturn) <- ntrials
  
  return(toreturn)
  
}





# 
# # T2 <- k_unfold(T_true,2)@data
# # Utrue2 <- irlba(T2,3,3)$u
# # T3 <- k_unfold(T_true,3)@data
# # Utrue3 <- irlba(T3,3,3)$u
# 
# #generate noise
# 
# 
# 
# 
# 
# #initialize true tensor:
# p <- 300
# #generate MMSBM
# 
# 
# 
# 
# 
# n1 = .3*p
# n2 = .2*p
# n3 = .5*p
# K= 3
# 
# W <- matrix(
#   c(1/2,1/2,-sqrt(2)/2
#     ,1/2,1/2,sqrt(2)/2
#     ,sqrt(2)/2,-sqrt(2)/2,0),
#   3,3
# )
# 
# B1 <- t(W) %*% diag(c(1.9,.1,.3),3,3) %*% W
# B2 <- t(W) %*% diag(c(1.9,-.1,-.3),3,3) %*% W 
# 
# Z <- matrix(c(
#   c(rep(1,n1),rep(0,n2),rep(0,n3)),
#   c(rep(0,n1),rep(1,n2),rep(0,n3)),
#   c(rep(0,n1),rep(0,n2),rep(1,n3))
# ),p,3
# 
# )
# L <- p/2
# P1 <- Z %*% B1 %*% t(Z)
# P2 <- Z %*% B2 %*% t(Z)
# sigma = 1
# 
# nMC <- 100
# dat_uhats <- matrix(0,nMC,K)
# dat_T <- matrix(0,nMC)
# for (mc in c(1:nMC)) {
#   print(paste(mc, "of",nMC))
#   T_hat <- array(0,dim = c(p,p,L))
#   T_true <- array(0,dim=c(p,p,L))
#   for (i in c(1:(L/2))) {
#     T_true[,,i] <- P1
#     T_hat[,,i] <- P1 + matrix(rnorm(p^2,0,sigma),p,p)
#   }
#   for (i in c((L/2 + 1):L)) {
#     T_true[,,i] <- P2
#     T_hat[,,i] <- P2 + matrix(rnorm(p^2,0,sigma),p,p)
#   }
#   T_hat <- as.tensor(T_hat)
#   T_true <- as.tensor(T_true)
#   
#   uhats <- HOOI_dd(T_hat,r = c(K,K,2))
#   #utrue <- HOOI_vanilla(T_true,r=c(K,K,2))
#   
#   T_est <- ttm(
#                 ttm(
#                   ttm(T_hat,uhats[[1]] %*% t(uhats[[1]]),1)
#                   ,uhats[[2]]%*%t(uhats[[2]]),2)
#                   ,uhats[[3]]%*%t(uhats[[3]]),3)
#   
#   
#   T1_svd <-irlba(k_unfold(T_true,1)@data,K)
#   T2_svd <-irlba(k_unfold(T_true,2)@data,K)
#   T3_svd <- irlba(k_unfold(T_true,3)@data,2)
#   
#   V1 <- T1_svd$v
#   V2 <- T2_svd$v
#   V3 <- T3_svd$v
#   
#   i1 <- 1
#   i2 <- 1
#   i3 <- 1
#   
#   variance_est <- sigma^2 * ( sum(abs(V1[i1,]^2)) + sum(abs(V2[i2,]^2)) + sum(abs(V3[i3,]^2)) )
#   
#   
#   dat_T[mc] <- T_est[1,1,1]@data - T_true[1,1,1]@data
#   dat_T[mc] <- dat_T[mc]/sqrt(variance_est)
#   
#   # 
#    utrue1 <- T1_svd$u
#    uhat1 <- uhats[[1]]
#    uhat1_proc <- uhat1 %*% procrustes(uhat1,utrue1)
#   # 
#    Lambda1 <- T1_svd$d
#    Lambda1 <- Lambda1/sigma
#   # 
#   dat_uhats[mc,] <- ((uhat1_proc - utrue1) %*% diag(Lambda1))[1,]
#   
# }
# 
# 
# save(list=ls(),file="4-13_sims.Rdata")
# load("4-13_sims.Rdata")
# dat_T <- as.data.frame(dat_T)
# gg <- ggplot(dat_T,aes(x = dat_T))
# gg <- gg + geom_histogram(bins=20, colour="black", 
#                           aes(y=..density..,fill=..count..))
# gg <- gg+guides(fill=FALSE)
# gg <- gg + scale_fill_gradient("Count", low="#DCDCDC", high="#7C7C7C")
# gg <- gg + stat_function(fun=dnorm,
#                          aes(linetype = "estimated"),color= "blue",
#                          args=list(mean=mean(dat$dat), 
#                                    sd=sd(dat$dat)))
# 
# gg <- gg + stat_function(aes(linetype = "theoretical"),color= "blue",fun=dnorm,
#                          args=list(mean=0, 
#                                    sd=1)) +
#   scale_linetype_manual("Overlaid Gaussian Density", values = c("dashed","solid"))
# gg
