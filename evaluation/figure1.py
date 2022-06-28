import nibabel as nib
import numpy as np
import ipdb
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcol
from matplotlib.colors import LinearSegmentedColormap



# ===== set colormap: `cm1` and `cm2`
# cm1 # black (bottom) - red - green (middle) - blue - white (top)
colors = [(0, 0, 0), (0.8, 0, 0.4), (0.6, 0.4, 0.2),
          (0.3, 0.8, 0), (0.3, 0.8, 0), (0.3, 0.8, 0), (0.3, 0.8, 0),
          (0.3, 1, 1), (1, 1, 1), (1, 1, 1), (1, 1, 1), (1, 1, 1)]
cm1 = LinearSegmentedColormap.from_list('my_list', colors, N=100)
# cm2 # black (bottom) - red - green (middle) - blue - white (top)
colors = [(1, 1, 1), (1, 1, 1), (1, 1, 1), (1, 1, 1), (0.3, 1, 1),
          (0.3, 0.8, 0), (0.3, 0.8, 0), (0.3, 0.8, 0), (0.3, 0.8, 0),
          (0.6, 0.4, 0.2), (0.8, 0, 0.4), (0, 0, 0)]
cm2 = LinearSegmentedColormap.from_list('my_list', colors, N=100)



# ====== load true values
# alpha true value
alphaT = np.loadtxt("./input/true/alpha_true.txt", dtype='f', delimiter=' ')
# beta true value
betaT = np.loadtxt("./input/true/beta_true.txt", dtype='f', delimiter=' ')



# ====== load estimates
# pick a replication (from 1 to 100)
rep = 5
# ----- BI-GMRF
# alpha
alpha_BG = np.loadtxt("./input/BI_GMRF/rep_alpha.txt", dtype='f')
alpha_BG = alpha_BG[rep,]
alpha_BG_mat = alpha_BG.reshape(64,64)
alpha_BG_mat = alpha_BG_mat.transpose()
# beta
beta_BG = np.loadtxt("./input/BI_GMRF/rep_beta.txt", dtype='f')
beta_BG = beta_BG[rep,]
beta_BG_mat = beta_BG.reshape(64,64)
beta_BG_mat = beta_BG_mat.transpose()
# ----- L2RM and TReg
# alpha
alpha_L2RM = np.loadtxt("./input/L2RM/rep_alpha.txt", dtype='f')
alpha_L2RM = alpha_L2RM[rep,]
alpha_L2RM_mat = alpha_L2RM.reshape(64,64)
alpha_L2RM_mat = alpha_L2RM_mat.transpose()
# beta rank=1
beta_TReg1 = np.loadtxt("./input/TReg/rep_beta_rk1.txt", dtype='f')
beta_TReg1 = beta_TReg1[rep,]
beta_TReg1_mat = beta_TReg1.reshape(64,64)
beta_TReg1_mat = beta_TReg1_mat.transpose()
# beta rank=2
beta_TReg2 = np.loadtxt("./input/TReg/rep_beta_rk2.txt", dtype='f')
beta_TReg2 = beta_TReg2[rep,]
beta_TReg2_mat = beta_TReg2.reshape(64,64)
beta_TReg2_mat = beta_TReg2_mat.transpose()
# beta rank=3
beta_TReg3 = np.loadtxt("./input/TReg/rep_beta_rk3.txt", dtype='f')
beta_TReg3 = beta_TReg3[rep,]
beta_TReg3_mat = beta_TReg3.reshape(64,64)
beta_TReg3_mat = beta_TReg3_mat.transpose()



# ====== plot figure 4
fig, axs = plt.subplots(3, 5, figsize=(8, 4.5), constrained_layout=True)
# min and max
minAl = -0.6
maxAl = 0.2
minBe = -0.2
maxBe = 0.6
minAB = -0.15
maxAB = 0.05
# true value
imAl = axs[0, 0].imshow(alphaT, vmax=maxAl, vmin=minAl, cmap=cm1)
imBe = axs[1, 0].imshow(betaT, vmax=maxBe, vmin=minBe, cmap=cm2)
imAB = axs[2, 0].imshow(alphaT*betaT, vmax=maxAB, vmin=minAB, cmap=cm1)
# BI-GMRF estimates
axs[0, 1].imshow(alpha_BG_mat, vmax=maxAl, vmin=minAl, cmap=cm1)
axs[1, 1].imshow(beta_BG_mat, vmax=maxBe, vmin=minBe, cmap=cm2)
axs[2, 1].imshow(alpha_BG_mat*beta_BG_mat, vmax=maxAB, vmin=minAB, cmap=cm1)
# L2RM and TReg rank=1 estimates
axs[0, 2].imshow(alpha_L2RM_mat, vmax=maxAl, vmin=minAl, cmap=cm1)
axs[1, 2].imshow(beta_TReg1_mat, vmax=maxBe, vmin=minBe, cmap=cm2)
axs[2, 2].imshow(alpha_L2RM_mat*beta_TReg1_mat, vmax=maxAB, vmin=minAB, cmap=cm1)
# L2RM and TReg rank=2 estimates
axs[0, 3].imshow(alpha_L2RM_mat, vmax=maxAl, vmin=minAl, cmap=cm1)
axs[1, 3].imshow(beta_TReg2_mat, vmax=maxBe, vmin=minBe, cmap=cm2)
axs[2, 3].imshow(alpha_L2RM_mat*beta_TReg2_mat, vmax=maxAB, vmin=minAB, cmap=cm1)
# L2RM and TReg rank=3 estimates
axs[0, 4].imshow(alpha_L2RM_mat, vmax=maxAl, vmin=minAl, cmap=cm1)
axs[1, 4].imshow(beta_TReg3_mat, vmax=maxBe, vmin=minBe, cmap=cm2)
axs[2, 4].imshow(alpha_L2RM_mat*beta_TReg3_mat, vmax=maxAB, vmin=minAB, cmap=cm1)
# add color bar
fig.colorbar(imAl, ax=[axs[0, -1]], shrink=0.8)
fig.colorbar(imBe, ax=[axs[1, -1]], shrink=0.8)
fig.colorbar(imAB, ax=[axs[2, -1]], shrink=0.8)

# change frame color and remove ticks for all
for i in range(3):
    for j in range(5):
        axs[i,j].spines['bottom'].set_color((0.98,0.98,0.98,0))
        axs[i,j].spines['top'].set_color((0.98,0.98,0.98,0))
        axs[i,j].spines['left'].set_color((0.98,0.98,0.98,0))
        axs[i,j].spines['right'].set_color((0.98,0.98,0.98,0))
        plt.setp(axs[i,j].get_xticklabels(), visible=False)
        plt.setp(axs[i,j].get_yticklabels(), visible=False)
        axs[i,j].tick_params(axis='both', which='both', length=0)

# background color
fig.patch.set_facecolor((0.98,0.98,0.98,0))

# save as .png
fig.savefig('./output/figure4.png', transparent=True)

# ipdb.set_trace()


