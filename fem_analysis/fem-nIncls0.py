from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
# import random
# from array import *
# import math
import numpy as np
# import os        # Operating system
# import shutil    # copying or moving files

####################################
## GLOBAL VARAIBLES AND ARGUMENTS ##
####################################
# constants (don't need to change this)
sName = '__profile__'   # sketch
mName = 'nIncls0'         # model
pName = 'Part-1'        # part
iName = 'Part-1-1'      # instance
foName = 'F-Output-1'   # field output

stepName = 'Step-1'
jName = mName

# variables (available for modifications)
model_width = 20
model_height = 10

dims = np.array([model_width, model_height, model_width])
sheetSize = 1.1*model_width if model_width>model_height else 1.1*model_height

diffCoeffMatrix = 1

########################################################################
# CREATE NEW MODEL OBJECT
########################################################################
mdb.Model(modelType=STANDARD_EXPLICIT, name=mName)
model = mdb.models[mName]

model.ConstrainedSketch(name=sName, sheetSize=sheetSize)
model.sketches[sName].rectangle(point1=(0, 0),
                                point2=(dims[0], dims[1]))
model.Part(dimensionality=THREE_D,
           name=pName,
           type=DEFORMABLE_BODY)
model.parts[pName].BaseSolidExtrude(depth=dims[2],
                                    sketch=model.sketches[sName])
part = model.parts[pName]

del model.sketches[sName]

########################################################################
# DEFINE MATERIAL PROPERTIES
########################################################################
# define material properties of matrix
model.Material(name='Matrix')
materMatrix = model.materials['Matrix']
materMatrix.Diffusivity(table=((diffCoeffMatrix, 0), ))
materMatrix.Solubility(table=((1, ), ))

########################################################################
# CREATE SECTIONS
########################################################################
model.HomogeneousSolidSection(material='Matrix',
                              name='Matrix',
                              thickness=None)

########################################################################
# ASSIGN SECTIONS TO PARTS
########################################################################
regionMatrix = Region(cells=part.cells.findAt(((.01, .01, .01),),))
part.SectionAssignment(offset=0.,
                       offsetField='',
                       offsetType=MIDDLE_SURFACE,
                       region=regionMatrix,
                       sectionName='Matrix',
                       thicknessAssignment=FROM_SECTION)

########################################################################
# CREATE INSTANCE (ASSEMBLY)
########################################################################
model.rootAssembly.DatumCsysByDefault(CARTESIAN)
model.rootAssembly.Instance(dependent=ON,
                            name=iName,
                            part=part)

########################################################################
# CREATE STEP
########################################################################
# model.MassDiffusionStep(initialInc=0.1,
#                         maxNumInc=int(1e8),
#                         name=stepName,
#                         previous='Initial',
#                         timeIncrementationMethod=FIXED,
#                         timePeriod=1000)

model.MassDiffusionStep(dcmax=100, # max allowable conc. change
                        initialInc=1e-3, # initial increment size
                        maxInc=10, # maximum increment size
                        minInc=1e-6, # minimum increment size
                        maxNumInc=int(1e8), # maximum number of increments
                        name=stepName,
                        previous='Initial',
                        timePeriod=200)
                        # end=0.1, # stop criteria (conc. change rate)

########################################################################
# SETUP BOUNDARY CONDITIONS
########################################################################
faceContSurf = model.rootAssembly.instances[iName].faces.findAt(((dims[0]/2,
                                                                  0,
                                                                  dims[2]/2),))
regionContSurf = Region(faces=faceContSurf)
model.ConcentrationBC(amplitude=UNSET,
                      createStepName=stepName,
                      distributionType=UNIFORM,
                      fieldName='',
                      fixed=OFF,
                      magnitude=1,
                      name='BC-1',
                      region=regionContSurf)

########################################################################
# SET ELEMENT TYPE
########################################################################
elemType1 = ElemType(elemCode=DC3D8, elemLibrary=STANDARD)
elemType2 = ElemType(elemCode=DC3D6, elemLibrary=STANDARD)
elemType3 = ElemType(elemCode=DC3D4, elemLibrary=STANDARD)

# elemType1 = ElemType(elemCode=DC3D20, elemLibrary=STANDARD)
# elemType2 = ElemType(elemCode=DC3D15, elemLibrary=STANDARD)
# elemType3 = ElemType(elemCode=DC3D10, elemLibrary=STANDARD)

part.setElementType(elemTypes=(elemType1, elemType2, elemType3),
                    regions=(part.cells.findAt(((.01, .01, .01), )), ))

########################################################################
# SEED MODEL
########################################################################
# global seeds
part.seedPart(constraint=FREE,
              deviationFactor=.1,
              minSizeFactor=.1,
              size=0.5)

########################################################################
# MESH MODEL
########################################################################
part.setMeshControls(elemShape=TET,
                     regions=part.cells.findAt(((.01, .01, .01), )),
                     technique=FREE)
part.generateMesh()

model.rootAssembly.regenerate()

########################################################################
# MODIFY OUTPUT REQUESTS
########################################################################
model.fieldOutputRequests[foName].setValues(position=AVERAGED_AT_NODES,
                                            variables=PRESELECT,
                                            timeInterval=1)

########################################################################
# CREATE JOB
########################################################################
mdb.Job(atTime=None,
        contactPrint=OFF,
        description='',
        echoPrint=OFF,
        explicitPrecision=SINGLE,
        getMemoryFromAnalysis=True,
        historyPrint=OFF,
        memory=90,
        memoryUnits=PERCENTAGE,
        model=mName,
        modelPrint=OFF,
        multiprocessingMode=DEFAULT,
        name=jName,
        nodalOutputPrecision=SINGLE,
        numCpus=1,
        numGPUs=0,
        queue=None,
        resultsFormat=ODB,
        scratch='',
        type=ANALYSIS,
        userSubroutine='',
        waitHours=0,
        waitMinutes=0)

# mdb.jobs[jName].submit(consistencyChecking=ON)
