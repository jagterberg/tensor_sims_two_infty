source("misc.R")
source("tensor_two_infty.R")
if (!require(pushoverr)) {
  install.packages("pushoverr")
  library(pushoverr)
}

print("got here!")
r <- 3
C <- array(runif(r^3),dim=rep(r,3))
C <- as.tensor(C)
p <- 300
ntrials <- 300
sigmas <- seq(5,50,1)

final_res <- two_infty_sim(p,r,C,sigmas,ntrials)
save(final_res,file = "sim_7-19.Rdata")
print("sims done.")
set_pushover_user(user = "utzcfs3wkwvzxi2h88w3dq1vixadsk")
set_pushover_app(token = "aznj99sdeso5cndn1mtrsb6ji1te9i")
pushover(message = "simulations finished!")
rowMeans(final_res)
