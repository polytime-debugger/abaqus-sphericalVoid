import numpy as np
import matplotlib.pyplot as plt

saveDir = '../'

########################################################################
# CONSTANTS
########################################################################
nTerms = 100 # approximate solution using num terms of sine

model_width = 20
model_height = 10
volBlock = model_width**2 * model_height

radiusIncl = 1
volIncl = 4/3 * radiusIncl**3 * np.pi

nBins = 50
binThickness = model_height/nBins
binLoc = np.linspace(0,model_height,nBins,endpoint=False) + binThickness/2

diffCoeff = 1
initialConc = 1

nFrames = 200
frameInterval = 1

nIncls_ls = np.array([0,80,160,240,300])

def effectDiffCoeff(volFrac, diffCoeff=diffCoeff):
    if volFrac>(1-1e-8):
        return diffCoeff
    else:
        return volFrac * diffCoeff / (1 - .5*np.log(volFrac))

########################################################################
# COMPUTE CONCENTRATION
########################################################################
concs = np.empty((len(nIncls_ls), nBins, nFrames))

for mIdx, nIncls in enumerate(nIncls_ls):
    volFrac = 1 - nIncls * volIncl / volBlock
    diffEff = effectDiffCoeff(volFrac)

    n_ls = np.arange(nTerms)
    a = -4 * initialConc / np.pi
    for frameIdx in range(nFrames):
        frameNum = frameIdx + 1
        t = frameNum * frameInterval
        for binIdx, x in enumerate(binLoc):
            X = a/(2*n_ls[1:]-1)*np.sin((n_ls[1:]-0.5)*np.pi*x/model_height)
            T = np.exp((-(n_ls[1:]-0.5)**2*np.pi**2*diffEff*t)/(model_height**2))
            concs[mIdx, binIdx, frameIdx] = initialConc + np.sum(X*T)

    ####################################################################
    # PLOT PROFILE
    ####################################################################
    for frameNum in 2**np.arange(1,8):
        frameIdx = frameNum - 1
        plt.plot(binLoc, concs[mIdx,:,frameIdx],
                 label='t={:d}'.format(frameNum*frameInterval))
    plt.title('concentration profile (analytical)\n'
              'volume fraction: {:.4f}'.format(volFrac))
    plt.xlabel('distance from source')
    plt.ylabel('concentration')
    plt.xlim(0,model_height)
    plt.ylim(0,initialConc)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    # plt.show()
    plt.savefig(saveDir+'profile-nIncls{}-analytical.png'.format(nIncls),
                dpi=300, bbox_inches='tight')
    plt.clf()

    ####################################################################
    # PLOT CHANGE
    ####################################################################
    # change at the center of last bin, not y=model_height
    t_ls = (np.arange(nFrames)+1) * frameInterval

    plt.plot(t_ls, concs[mIdx,-1])
    plt.title('change of concentration at y=9.9 over time (analytical)\n'
              'volume fraction: {:.4f}'.format(volFrac))
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.xlim(0,t_ls[-1])
    plt.ylim(0,initialConc)
    plt.grid()
    plt.tight_layout()
    # plt.show()
    plt.savefig(saveDir+'change-nIncls{}-analytical.png'.format(nIncls),
                dpi=300, bbox_inches='tight')
    plt.clf()

########################################################################
# PLOT CHANGE (ALL)
########################################################################
t_ls = (np.arange(nFrames)+1) * frameInterval
for mIdx, nIncls in enumerate(nIncls_ls):
    plt.plot(t_ls, concs[mIdx,-1], label=nIncls)
plt.title('change of concentration at y=9.9 over time\n(analytical)')
plt.xlabel('time')
plt.ylabel('concentration')
plt.xlim(0,t_ls[-1])
plt.ylim(0,initialConc)
plt.legend(title='volume fraction')
plt.grid()
plt.tight_layout()
plt.savefig(saveDir+'change-all-analytical.png', dpi=300, bbox_inches='tight')


########################################################################
# SAVE DATA TO DISK
########################################################################
np.save(saveDir+'concs-analytical-5MODELSx50BINSx200FRAMES', concs)
