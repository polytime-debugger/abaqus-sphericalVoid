import numpy as np
import csv
import matplotlib.pyplot as plt

saveDir = '../'
dataDir = '../../fem/data/'
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

nIncls = 240

mName_ls = ['nIncls240m1','nIncls240','nIncls240m2']
meshSize_ls =  [0.5, 0.4, 0.3]

####################################################################
# GET CONCENTRATION (FROM FILES)
####################################################################
concs = np.empty((len(mName_ls), nBins, nFrames))

for mIdx, mName in enumerate(mName_ls):
    print('working on model "{}"\n'.format(mName))

    coordsFile = mName+'-coords.csv'
    binNodeFile = mName+'-nodes-bin50.csv'
    concFile = mName+'-concs.npy'
    
    volFrac = 1 - nIncls * volIncl / volBlock
    
    ####################################################################
    # READ DATA
    ####################################################################
    # coordinates
    data = np.genfromtxt(dataDir+coordsFile, delimiter=' ')
    indices = data[:,[0]].astype(np.int32)
    coords = data[:,1:]
    print('indices.shape:', indices.shape)
    print('coords.shape:', coords.shape)

    # concentration at nodes
    concsNodes = np.load(dataDir+concFile)
    print('concsNodes.shape:', concsNodes.shape, end='\n\n')

    # node indices in bins
    binNodeNums_ls = []
    with open(dataDir+binNodeFile, 'r') as fin:
        reader = csv.reader(fin)
        for row in reader:
            binNodeNums = np.array([int(idx) for idx in row], dtype=np.int32)
            binNodeNums_ls.append(binNodeNums)

    nNodes = indices.size
    print('number of nodes:', nNodes)
    nFrames = concsNodes.shape[1]-1
    print('total frames: {} (={}-1)'.format(nFrames, concsNodes.shape[1]))
    nBins = len(binNodeNums_ls)
    print('number of bins:', nBins, end='\n\n')
    
    ####################################################################
    # COMPUTE AVERAGE OVER ALL NODES IN BINS
    ####################################################################
    for binIdx, binNodeNums in enumerate(binNodeNums_ls):
        concs[mIdx,binIdx] = np.average(concsNodes[binNodeNums-1,1:], axis=0)

    ####################################################################
    # PLOT PROFILE
    ####################################################################
    for frameNum in 2**np.arange(1,8):
        frameIdx = frameNum - 1
        plt.plot(binLoc, concs[mIdx,:,frameIdx],
                 label='t={:d}'.format(frameNum*frameInterval))
    plt.title('concentration profile (FEM)\n'
              'volume fraction: {:.4f}, mesh size: {:.2f}'
              ''.format(volFrac, meshSize_ls[mIdx]))
    plt.xlabel('distance from source')
    plt.ylabel('concentration')
    plt.xlim(0,model_height)
    plt.ylim(0,initialConc)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    # plt.show()
    plt.savefig(saveDir+'profile-{}-fem240.png'.format(mName),
                dpi=300, bbox_inches='tight')
    plt.clf()

    ####################################################################
    # PLOT CHANGE
    ####################################################################
    # change at the center of last bin, not y=model_height
    t_ls = (np.arange(nFrames)+1) * frameInterval

    plt.plot(t_ls, concs[mIdx,-1])
    plt.title('change of concentration at y=9.9 over time (FEM)\n'
              'volume fraction: {:.4f}, mesh size: {:.2f}'
              ''.format(volFrac, meshSize_ls[mIdx]))
    plt.xlabel('time')
    plt.ylabel('concentration')
    plt.xlim(0,t_ls[-1])
    plt.ylim(0,initialConc)
    plt.grid()
    plt.tight_layout()
    # plt.show()
    plt.savefig(saveDir+'change-{}-fem240.png'.format(mName),
                dpi=300, bbox_inches='tight')
    plt.clf()

########################################################################
# PLOT CHANGE (ALL)
########################################################################
t_ls = (np.arange(nFrames)+1) * frameInterval
for mIdx, mName in enumerate(mName_ls):
    plt.plot(t_ls, concs[mIdx,-1], label=meshSize_ls[mIdx])
plt.title('change of concentration at y=9.9 over time\n'
          'convergence study')
plt.xlabel('time')
plt.ylabel('concentration')
plt.xlim(0,t_ls[-1])
plt.ylim(0,initialConc)
plt.legend(title='mesh size')
plt.grid()
plt.tight_layout()
# plt.show()
plt.savefig(saveDir+'change-all-fem240.png', dpi=300, bbox_inches='tight')


########################################################################
# SAVE DATA TO DISK
########################################################################
np.save(saveDir+'concs-fem240-3MESHESx50BINSx196FRAMES', concs)
