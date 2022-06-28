### create Table 1 ###
# RMSEs of imaging-related coefficients

library(xtable)
library(stargazer)
options(scipen = 100) # remove scientific notation
nDec=4 # number of decimal places
V = 64*64 # total number of pixels
nRep = 100 # number of replications

# function for summarize the result of one scenario
# scenario = c("square", "cross", "circle", "triangle", "butterfly")
summ = function(scenario)
{
  path = paste("input/", scenario, "/", sep="")
  # read true values
  alphaT = as.vector(as.matrix(read.table(paste(path, "true/alpha_true.txt", sep=""), header = F)))
  betaT = as.vector(as.matrix(read.table(paste(path, "true/beta_true.txt", sep=""), header = F)))
  abT = alphaT*betaT
  # read BI-GMRF estimates and compute RMSE
  alpha_BG = as.matrix(read.table(paste(path, "BI_GMRF/rep_alpha.txt", sep=""), header = F))
  alpha_BG_rmse = sqrt(sum((t(alpha_BG) - alphaT)^2)/(nRep*V))
  beta_BG = as.matrix(read.table(paste(path, "BI_GMRF/rep_beta.txt", sep=""), header = F))
  beta_BG_rmse = sqrt(sum((t(beta_BG) - betaT)^2)/(nRep*V))
  ab_BG = alpha_BG*beta_BG
  ab_BG_rmse = sqrt(sum((t(ab_BG) - abT)^2)/(nRep*V))
  # read L2RM estimates and compute RMSE
  alpha_L2RM = as.matrix(read.table(paste(path, "L2RM/rep_alpha.txt", sep=""), header = F))
  alpha_L2RM_rmse = sqrt(sum((t(alpha_L2RM) - alphaT)^2)/(nRep*V))
  # read TReg estimates (rank=1,2,3) and compute RMSE
  beta_TReg1 = as.matrix(read.table(paste(path, "TReg/rep_beta_rk1.txt", sep=""), header = F))
  beta_TReg2 = as.matrix(read.table(paste(path, "TReg/rep_beta_rk2.txt", sep=""), header = F))
  beta_TReg3 = as.matrix(read.table(paste(path, "TReg/rep_beta_rk3.txt", sep=""), header = F))
  beta_TReg1_rmse = sqrt(sum((t(beta_TReg1) - betaT)^2)/(nRep*V))
  beta_TReg2_rmse = sqrt(sum((t(beta_TReg2) - betaT)^2)/(nRep*V))
  beta_TReg3_rmse = sqrt(sum((t(beta_TReg3) - betaT)^2)/(nRep*V))
  # combine L2RM alpha and TReg beta to ab_LR and compute RMSE
  ab_LR1 = alpha_L2RM*beta_TReg1
  ab_LR2 = alpha_L2RM*beta_TReg2
  ab_LR3 = alpha_L2RM*beta_TReg3
  ab_LR1_rmse = sqrt(sum((t(ab_LR1) - abT)^2)/(nRep*V))
  ab_LR2_rmse = sqrt(sum((t(ab_LR2) - abT)^2)/(nRep*V))
  ab_LR3_rmse = sqrt(sum((t(ab_LR3) - abT)^2)/(nRep*V))
  # array RMSEs
  tab = array("", c(3, 7))
  tab[,1] = c("alpha","beta","ab")
  tab[,2] = scenario
  tab[,3] = round(c(alpha_BG_rmse, beta_BG_rmse, ab_BG_rmse), nDec)
  tab[1,4] = round(alpha_L2RM_rmse, nDec)
  tab[2:3,5] = round(c(beta_TReg1_rmse, ab_LR1_rmse), nDec)
  tab[2:3,6] = round(c(beta_TReg2_rmse, ab_LR2_rmse), nDec)
  tab[2:3,7] = round(c(beta_TReg3_rmse, ab_LR3_rmse), nDec)

  # return tab
  return(tab)
}

# RMSE of 5 scenarios
tab_square = summ("square")
tab_cross = summ("cross")
tab_circle = summ("circle")
tab_triangle = summ("triangle")
tab_butterfly = summ("butterfly")

# combine to table
table = rbind(tab_square[1,], tab_cross[1,], tab_circle[1,], tab_triangle[1,], tab_butterfly[1,], # alpha
              tab_square[2,], tab_cross[2,], tab_circle[2,], tab_triangle[2,], tab_butterfly[2,], # beta
              tab_square[3,], tab_cross[3,], tab_circle[3,], tab_triangle[3,], tab_butterfly[3,]) # ab

# save as latex table
table_df = as.data.frame(table)
colnames(table_df) = c("para", "scenarios", "BI-GMRF", "L2RM", "TReg rank=1", "TReg rank=2", "TReg rank=3")
xtable(table_df, type = "latex", sink("./output/table1.tex"))















