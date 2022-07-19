source("misc.R")
source("tensor_two_infty.R")

r <- 3
C <- array(runif(r^3),dim=rep(r,3))
C <- as.tensor(C)
p <- 100
ntrials <- 100
sigmas <- seq(1,50,5)

final_res <- tensor_two_infty(p,r,C,sigmas,ntrials)
save(final_res,file = "sim_7-19.Rdata")

