# from abaqus import *
from odbAccess import *
import numpy as np

########################################################################
# constants
########################################################################
pName = 'Part-1'
iName = 'Part-1-1'
sName = 'Step-1'

iNameODB = 'PART-1-1'

model_width = 20
model_height = 10

saveDir = '../data/'

nIncls = [0,80,160,240,320]

# mName_ls = ['nIncls{}'.format(n) for n in nIncls]
# mName_ls = ['nIncls240m1', 'nIncls240m2']
mName_ls = ['nIncls240m2']

for mName in mName_ls:

    print '\nSaving CONC of model "{}"'.format(mName)

    ########################################################################
    # GET NODES
    ########################################################################
    model = mdb.models[mName]
    nodes = model.rootAssembly.instances[iName].nodes

    xmin, ymin, zmin = 0-1, 0-1, 0-1
    xmax, ymax, zmax = model_width+1, model_height+1, model_width+1
    nodesAll = nodes.getByBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)

    nNodes = len(nodesAll)
    print nNodes, 'nodes in model'

    labelsAll = np.array([[node.label] for node in nodesAll])
    coordsAll = np.array([node.coordinates for node in nodesAll])

    ########################################################################
    # GET CONCENTRATION
    ########################################################################
    odb = openOdb(path='{}.odb'.format(mName))

    nFrames = len(odb.steps[sName].frames)
    print nFrames, 'frames in record'

    labels = np.empty((nNodes, 1))
    concs = np.empty((nNodes, nFrames))

    for frameIdx, frame in enumerate(odb.steps[sName].frames):
        if frameIdx%10 == 0:
            print 'Frame {}'.format(frameIdx)

        concFrame = frame.fieldOutputs['CONC']

        for idx in range(nNodes):

            labels[idx] = concFrame.values[idx].nodeLabel
            concs[idx, frameIdx] = concFrame.values[idx].data

    # dataFrame = np.hstack((labels, concs))
    # np.savetxt(saveDir+mName+'-CONC.csv', dataFrame,
    #            fmt=['%6d']+nFrames*['%.18e'],
    #            delimiter=' ')
    np.save(saveDir+mName+'-concs', concs)

    odb.close()
