from abaqus import *
import numpy as np
import csv

########################################################################
# constants
########################################################################
iName = 'Part-1-1'

model_width = 20
model_height = 10

nBins = 50

saveDir = '../data/'

nIncls = [0,80,160,240,320]

# mName_ls = ['nIncls{}'.format(n) for n in nIncls]
# mName_ls = ['nIncls240m1', 'nIncls240m2']
mName_ls = ['nIncls240m2']

for mName in mName_ls:

    print '\nSaving coordinates of model "{}"'.format(mName)

    ########################################################################
    # GET NODES
    ########################################################################
    model = mdb.models[mName]
    nodes = model.rootAssembly.instances[iName].nodes

    xmin, ymin, zmin = 0-1, 0-1, 0-1
    xmax, ymax, zmax = model_width+1, model_height+1, model_width+1
    nodesAll = nodes.getByBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)

    nodesBin_ls = []
    binThickness = float(model_height)/float(nBins)

    for i in range(nBins):
        xmin, ymin, zmin = 0, i*binThickness, 0
        xmax, ymax, zmax = model_width, (i+1)*binThickness, model_width
        nodesBin = nodes.getByBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)
        nodesBin_ls.append(nodesBin)

    ########################################################################
    # SAVE TO CSV
    ########################################################################
    nNodes = len(nodesAll)
    print 'total num of nodes:', nNodes


    # Save all node labels and coordinates
    labels = np.array([[node.label] for node in nodesAll])
    coords = np.array([node.coordinates for node in nodesAll])

    dataAll = np.hstack((labels,coords))
    np.savetxt(saveDir+mName+'-coords.csv', dataAll,
               fmt=['%6d','%.18e','%.18e','%.18e'],
               delimiter=' ')


    # save each bin of node indices as a line
    with open(saveDir+mName+'-nodes-bin{}.csv'.format(nBins), 'wb') as fout:
        writer = csv.writer(fout)

        count = 0
        for nodesBin in nodesBin_ls:
            node_ls = [node.label for node in nodesBin]
            count += len(node_ls)

            writer.writerow(node_ls)
