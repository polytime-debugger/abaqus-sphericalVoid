import numpy as np
import matplotlib.pyplot as plt
# import matplotlib as mpl
# mpl.rcParams['text.usetex'] = True

'''
Plot error against time. 
One graph for each value of volume fraction.
'''

saveDir = '../'
dataDir = '../'
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
    if volFrac<(1-1e-8):
        return diffCoeff
    else:
        return volFrac * diffCoeff / (1 - .5*np.log(volFrac))

########################################################################
# READ DATA
########################################################################
analytical = np.load(dataDir+'concs-analytical-5MODELSx50BINSx200FRAMES.npy')
numerical  = np.load(dataDir+'concs-fem-5MODELSx50BINSx200FRAMES.npy')

print('\nanalytical.shape:', analytical.shape)
print('numerical.shape:', numerical.shape)
nBins = analytical.shape[1]
nFrames = analytical.shape[2]
print('number of bins:', nBins)
print('number of frames:', nFrames)

for mIdx, nIncls in enumerate(nIncls_ls):
    mName = 'nIncls{}'.format(nIncls)
    print('working on model "{}"\n'.format(mName))

    volFrac = 1 - nIncls*volIncl/volBlock
    diffEff = effectDiffCoeff(volFrac)

    ####################################################################
    # PLOT CHANGE AT BIN #50 (ANALYTICAL VS FEM )
    ####################################################################
    t_ls = (np.arange(nFrames)+1) * frameInterval

    fig, ax1 = plt.subplots()
    ax1.plot(t_ls, analytical[mIdx,-1], label='analytical')
    ax1.plot(t_ls, numerical[mIdx,-1], label='FEM')
    ax1.set_xlabel('time')
    ax1.set_ylabel('concentration')
    ax1.set_xlim(0,t_ls[-1])
    ax1.set_ylim(0,initialConc)
    ax1.legend(loc='right')
    ax1.grid()
    ax1.tick_params(axis='y', colors='k')

    ax2 = ax1.twinx()
    ax2.plot(t_ls,
             100*np.abs(analytical[mIdx,-1]-numerical[mIdx,-1])/analytical[mIdx,-1],
             label='error', c='r')
    ax2.set_ylabel('error (%)')
    ax2.set_ylim(0,100)
    ax2.tick_params(axis='y', colors='r')
    ax2.yaxis.label.set_color('red')

    fig.suptitle('error (analytic vs FEM) in concentration at y=9.9 over time\n'
                 'volume fraction: {:.4f}, $D_{{effective}}$: {:.4f}'
                 ''.format(volFrac, diffEff))
    fig.tight_layout()
    # plt.show()
    plt.savefig(saveDir+'error-change-nIncls{}'.format(nIncls),
                dpi=300,
                bbox_inches='tight')
    plt.clf()

    ####################################################################
    # PLOT PROFILE (ANALYTICAL VS FEM )
    ####################################################################
    # error axis limited
    for frameNum in 2**np.arange(1,8):
        frameIdx = frameNum - 1
        anaFrame = analytical[mIdx,:,frameIdx]
        numFrame = numerical[mIdx,:,frameIdx]
        plt.plot(binLoc,
                 100*np.abs(anaFrame-numFrame)/anaFrame,
                 label='t={:d}'.format(frameNum*frameInterval))
    plt.title('error (analytic vs FEM) in concentration profile')
    plt.xlabel('distance from source')
    plt.ylabel('error (%)')
    plt.xlim(0,model_height)
    plt.ylim(0,100)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    # plt.show()
    plt.savefig(saveDir+'error-profile-nIncls{}-limited.png'.format(nIncls),
                dpi=300, bbox_inches='tight')
    plt.clf()

    # error axis unlimited
    for frameNum in 2**np.arange(1,8):
        frameIdx = frameNum - 1
        anaFrame = analytical[mIdx,:,frameIdx]
        numFrame = numerical[mIdx,:,frameIdx]
        plt.plot(binLoc,
                 100*np.abs(anaFrame-numFrame)/anaFrame,
                 label='t={:d}'.format(frameNum*frameInterval))
    plt.title('error (analytic vs FEM) in concentration profile')
    plt.xlabel('distance from source')
    plt.ylabel('error (%)')
    plt.xlim(0,model_height)
    # plt.ylim(0,100)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    # plt.show()
    plt.savefig(saveDir+'error-profile-nIncls{}.png'.format(nIncls),
                dpi=300, bbox_inches='tight')
    plt.clf()
    # break
