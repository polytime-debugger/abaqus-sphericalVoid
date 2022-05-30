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

for nIncls in [240]: #[80,160,240,320]:
    ####################################
    ## GLOBAL VARAIBLES AND ARGUMENTS ##
    ####################################
    # constants (don't need to change this)
    sName = '__profile__'               # sketch
    mName = 'nIncls{}m2'.format(nIncls)   # model
    pName = 'Part-1'                    # part
    iName = 'Part-1-1'                  # instance
    foName = 'F-Output-1'               # field output

    stepName = 'Step-1'
    jName = mName

    # variables (available for modifications)
    radiusIncl = 1
    model_width = 20
    model_height = 10

    dims = np.array([model_width, model_height, model_width])
    sheetSize = 1.1*model_width if model_width>model_height else 1.1*model_height

    diffCoeffMatrix = 1


    ########################################################################
    # GENERATE INCLUSION LOCATION
    ########################################################################
    # list of radii
    radii = np.ones((nIncls, )) * radiusIncl

    # list of coordinates
    coords = np.zeros((nIncls, 3))

    count = 0
    iterCnt = 0
    while (count < nIncls):

        newLoc = np.random.rand(3) * (dims-2*radiusIncl) + np.ones((3,))*radiusIncl

        notIntersecting = True
        for idx in range(count):

            checkDist = 2 * radiusIncl

            if np.linalg.norm(newLoc-coords[idx]) < checkDist:
                notIntersecting = False
                iterCnt += 1
                break

        if notIntersecting:
            # print count, radii[count]
            coords[count] = newLoc
            count += 1
            iterCnt = 0
        elif iterCnt > 1e5:
            print 'Cannot find location to insert inclusion after 10000 tries'
            radii = radii[:count]
            coords = coords[:count]
            nIncls = count
            break

    print 'nIncls:', count

    ########################################################################
    # GENERATE ARGS FOR LATER USE
    ########################################################################
    # list of inlusion region coordinates
    radiusMid = radii/2
    coordChange = radiusMid / (3**.5)

    coordsUpper = coords + coordChange[:, np.newaxis]
    coordsLower = coords - coordChange[:, np.newaxis]

    coordsIncl = tuple(map(tuple, np.vstack((coordsUpper, coordsLower))))
    coordsIncl = tuple((coord,) for coord in coordsIncl)

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
    # CREATE PARTITIONS FOR INCLUSIONS
    ########################################################################
    # datum id placeholders
    dpYZ_Idxs = np.zeros((nIncls, ), dtype=int) # object)
    dpXZ_Idxs = np.zeros((nIncls, ), dtype=int) # object)
    dpXY_Idxs = np.zeros((nIncls, ), dtype=int) # object)

    for inclIdx in range(radii.size):
        # print inclIdx
        x, y, z = coords[inclIdx]

        radius = radii[inclIdx]

        ####################################################################
        ## Create Datum Planes
        ####################################################################
        dpYZ = part.DatumPlaneByPrincipalPlane(offset=x, principalPlane=YZPLANE)
        dpXZ = part.DatumPlaneByPrincipalPlane(offset=y, principalPlane=XZPLANE)
        dpXY = part.DatumPlaneByPrincipalPlane(offset=z, principalPlane=XYPLANE)

        datumYZ = part.datums[dpYZ.id]
        datumXZ = part.datums[dpXZ.id]
        datumXY = part.datums[dpXY.id]

        dpYZ_Idxs[inclIdx] = dpYZ.id
        dpXZ_Idxs[inclIdx] = dpXZ.id
        dpXY_Idxs[inclIdx] = dpXY.id

        ####################################################################
        ## Cut out sphere
        ####################################################################
        datumPlane = datumYZ
        sketchUpEdge = part.edges.findAt((0, dims[1]/2, 0), )

        # 4x3 matrix representing the transformation from sketch to part coords
        transform = part.MakeSketchTransform(sketchPlane=datumPlane,
                                             sketchPlaneSide=SIDE1,
                                             sketchUpEdge=sketchUpEdge,
                                             sketchOrientation=RIGHT,
                                             origin=(x, y, z))
        model.ConstrainedSketch(gridSpacing=2,
                                name=sName,
                                sheetSize=sheetSize,
                                transform=transform)
        # switch workspace to the sketch plane (from 3D to 2D)
        part.projectReferencesOntoSketch(sketch=model.sketches[sName],
                                         filter=COPLANAR_EDGES)

        # sketch semicircle for revolution-deletion
        model.sketches[sName].ArcByCenterEnds(center=(0,0),
                                              direction=CLOCKWISE,
                                              point1=(0, radius),
                                              point2=(0, -radius))
        model.sketches[sName].Line(point1=(0, radius),
                                   point2=(0, -radius))

        model.sketches[sName].ConstructionLine(angle=90,
                                               point1=(0, 0))
        centerline = model.sketches[sName].geometry.findAt((0, 111114.052), )
        model.sketches[sName].assignCenterline(line=centerline)

        # cut
        part.CutRevolve(angle=360,
                        flipRevolveDirection=OFF,
                        sketch=model.sketches[sName],
                        sketchOrientation=RIGHT,
                        sketchPlane=datumPlane,
                        sketchPlaneSide=SIDE1,
                        sketchUpEdge=sketchUpEdge)

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
    faceContSurf = model.rootAssembly.instances[iName].faces.findAt(((dims[0]/2, 0, dims[2]/2), ))
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
                  size=0.3)

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

    mdb.jobs[jName].submit(consistencyChecking=ON)
