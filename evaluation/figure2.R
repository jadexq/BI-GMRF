#################################################################################################
# ROC curves for butterfly shape, low signal-to-noisy ratio scenario:
#             alpha: compare BI-GMRF, BI-GMRF-0, and L2RM
#             beta:  compare BI-GMRF, BI-GMRF-0, and TReg rank=1 (corresponding to the best BIC)
#################################################################################################
library(plotly)
library(multiplex)
library(processx)

# To save the plotly images to files, the export approach is used. For a more updated approach, orca can be used. 
library(webshot)
webshot::install_phantomjs()

options("scipen"=100, "digits"=5)



# ========== constants
V = 64*64 # number of pixels
N = 2000 # number of grids for ROC curve
t = list(size = 30) # font size in figures



# ========== function computing TPR and FPR
# input: vecEst - estimated imaging coefficient (vecterized)
#        vecT - true imaging coefficient (vecterized)
#        N - number of grids
# output: array = (TPR | FPR)
rates <- function(vecEst, vecT, N)
{
  thre1 = 0 # min
  thre2 = max(abs(vecEst)) # max
  gridROC = seq(thre1, thre2, length.out = N) # grid for ROC
  TPR = rep(NA, N) # true positive rate
  FPR = rep(NA, N) # false positive rate
  for (i in 1:N)
  {
    cut = gridROC[i]
    TP = sum(abs(vecT[abs(vecEst) > cut]) > 0)
    FN = sum(abs(vecT[abs(vecEst) < cut]) > 0)
    FP = sum(abs(vecT[abs(vecEst) > cut]) == 0)
    TN = sum(abs(vecT[abs(vecEst) < cut]) == 0)
    TPR[i] = TP/(TP+FN)
    FPR[i] = FP/(FP+TN)
  }
  rates = cbind(TPR, FPR)
  return(rates)
}



# ========== read true values
# alpha
alphaT = as.vector(as.matrix(read.table("./input/true/alpha_true.txt", header = F)))
# beta
betaT = as.vector(as.matrix(read.table("./input/true/beta_true.txt", header = F)))



# ========== read estimates
# BI-GMRF
alpha_BG = as.vector(as.matrix(read.table("./input/BI_GMRF/rep_alpha.txt", header = F)))
beta_BG = as.vector(as.matrix(read.table("./input/BI_GMRF/rep_beta.txt", header = F)))
# BI-GMRF-0
alpha_BG0 = as.vector(as.matrix(read.table("./input/BI_GMRF0/rep_alpha.txt", header = F)))
beta_BG0 = as.vector(as.matrix(read.table("./input/BI_GMRF0/rep_beta.txt", header = F)))
# L2RM
alpha_L2RM = as.vector(as.matrix(read.table("./input/L2RM/rep_alpha.txt", header = F)))
# TReg rank=1,2,3 (rank=1 has the best BIC, see BIC.jpg)
beta_TReg1 = as.vector(as.matrix(read.table("./input/TReg/rep_beta_rk1.txt", header = F)))
#beta_TReg2 = as.vector(as.matrix(read.table("./input/TReg/rep_beta_rk2.txt", header = F)))
#beta_TReg3 = as.vector(as.matrix(read.table("./input/TReg/rep_beta_rk3.txt", header = F)))



# ========== ROC curves
# ----- alpha
# compute TPR and FPR
rate_alpha_BG = rates(alpha_BG, alphaT, N) # BI-GMRF
rate_alpha_BG0 = rates(alpha_BG0, alphaT, N) # BI-GMRF-0
rate_alpha_L2RM = rates(alpha_L2RM, alphaT, N) # L2RM
# plot
figA <- plot_ly(x = ~rate_alpha_L2RM[,2], y = ~rate_alpha_L2RM[,1], type = 'scatter', mode = 'lines', name = "L2RM", line = list(color = 'rgb(205, 0, 0)'))
figA <- figA %>% add_trace(x = ~rate_alpha_BG0[,2], y = ~rate_alpha_BG0[,1], type = 'scatter', mode = 'lines', name = "BI-GMRF-0", line = list(color = 'rgb(132, 112, 255)', dash = 'dot'))
figA <- figA %>% add_trace(x = ~rate_alpha_BG[,2], y = ~rate_alpha_BG[,1], type = 'scatter', mode = 'lines', name = "BI-GMRF", line = list(color = 'rgb(0, 205, 0)'))
figA <- figA %>% layout(title = '',
                        font=t,
                        xaxis = list(title = '', showgrid = F),
                        yaxis = list (title = '', showgrid = F),
                        paper_bgcolor='rgba(0, 0, 0, 0)',
                        plot_bgcolor='rgba(0, 0, 0, 0)',
                        width=500, 
                        height=500,
                        legend = list(x = 0.5, y = 0.1)) 
# save as .png
#orca(figA, file="./output/figA.png")
export(figA, file="./output/figA.png")

# ----- beta
# compute TPR and FPR
rate_beta_BG = rates(beta_BG, betaT, N) # BI-GMRF
rate_beta_BG0 = rates(beta_BG0, betaT, N) # BI-GMRF-0
rate_beta_TReg1 = rates(beta_TReg1, betaT, N) # TReg rank=1
# plot
figB <- plot_ly(x = ~rate_beta_TReg1[,2], y = ~rate_beta_TReg1[,1], type = 'scatter', mode = 'lines', name = "TReg rank=1", line = list(color = 'rgb(240,128,128)'))
figB <- figB %>% add_trace(x = ~rate_beta_BG0[,2], y = ~rate_beta_BG0[,1], type = 'scatter', mode = 'lines', name = "BI-GMRF-0", line = list(color = 'rgb(132, 112, 255)', dash = 'dot'))
figB <- figB %>% add_trace(x = ~rate_beta_BG[,2], y = ~rate_beta_BG[,1], type = 'scatter', mode = 'lines', name = "BI-GMRF", line = list(color = 'rgb(0,205,0)'))
figB <- figB %>% layout(title = '',
                        font=t,
                        xaxis = list(title = '', showgrid = F),
                        yaxis = list (title = '', showgrid = F),
                        paper_bgcolor='rgba(0, 0, 0, 0)',
                        plot_bgcolor='rgba(0, 0, 0, 0)',
                        width=500,
                        height=500,
                        legend = list(x = 0.42, y = 0.1))
# save as .png
#orca(figB, file="./output/figB.png")
export(figB, file="./output/figB.png")

# ----- alpha*beta
# compute TPR and FPR
rate_ab_BG = rates(alpha_BG*beta_BG, alphaT*betaT, N) # BI-GMRF
rate_ab_BG0 = rates(alpha_BG0*beta_BG0, alphaT*betaT, N) # BI-GMRF-0
rate_ab_LR1 = rates(alpha_L2RM*beta_TReg1, alphaT*betaT, N) # L2RM + TReg rank=1
# plot
figAB <- plot_ly(x = ~rate_ab_LR1[,2], y = ~rate_ab_LR1[,1], type = 'scatter', mode = 'lines', name = "L2RM+TReg rank=1", line = list(color = 'rgb(240,128,128)'))
figAB <- figAB %>% add_trace(x = ~rate_ab_BG0[,2], y = ~rate_ab_BG0[,1], type = 'scatter', mode = 'lines', name = "BI-GMRF-0", line = list(color = 'rgb(132, 112, 255)', dash = 'dot'))
figAB <- figAB %>% add_trace(x = ~rate_ab_BG[,2], y = ~rate_ab_BG[,1], type = 'scatter', mode = 'lines', name = "BI-GMRF", line = list(color = 'rgb(0,205,0)'))
figAB <- figAB %>% layout(title = '',
                          font=t,
                          xaxis = list(title = '', showgrid = F),
                          yaxis = list (title = '', showgrid = F),
                          paper_bgcolor='rgba(0, 0, 0, 0)',
                          plot_bgcolor='rgba(0, 0, 0, 0)',
                          width=500,
                          height=500,
                          legend = list(x = 0.05, y = 0.1))
# save as .png
#orca(figAB, file="./output/figAB.png")
export(figAB, file="./output/figAB.png")


