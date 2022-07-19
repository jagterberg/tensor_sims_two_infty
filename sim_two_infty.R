source("misc.R")
source("tensor_two_infty.R")

print("got here!")
r <- 3
C <- array(runif(r^3),dim=rep(r,3))
C <- as.tensor(C)
p <- 100
ntrials <- 100
sigmas <- seq(1,50,5)

final_res <- two_infty_sim(p,r,C,sigmas,ntrials)
save(final_res,file = "sim_7-19.Rdata")
print("sims done.")
