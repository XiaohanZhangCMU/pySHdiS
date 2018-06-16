from abaqus import *
from abaqusConstants import *
#session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=338.666656494141, 
#    height=154.820419311523)
#session.viewports['Viewport: 1'].makeCurrent()
#session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
#session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
#    referenceRepresentation=ON)
Mdb()
#: A new model database has been created.
#: The model "Model-1" has been created.
#session.viewports['Viewport: 1'].setValues(displayedObject=None)
mdb.models.changeKey(fromName='Model-1', toName='Beam')

#session.viewports['Viewport: 1'].setValues(displayedObject=None)
s = mdb.models['Beam'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(-50, 50), point2=(50.0, -50))
p = mdb.models['Beam'].Part(name='myBeam', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Beam'].parts['myBeam']
p.BaseSolidExtrude(sketch=s, depth=100.0)
s.unsetPrimaryObject()
p = mdb.models['Beam'].parts['myBeam']
#session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Beam'].sketches['__profile__']
#session.viewports['Viewport: 1'].view.setValues(nearPlane=86.8842, 
#    farPlane=164.567, width=83.5614, height=38.0978, cameraPosition=(111.922, 
#    46.289, -12.7409), cameraUpVector=(-0.694152, 0.712221, -0.104376), 
#    cameraTarget=(-0.915448, -1.98662, 9.77707))
#session.viewports['Viewport: 1'].view.setValues(nearPlane=88.3632, 
#    farPlane=163.269, width=84.9839, height=38.7463, cameraPosition=(94.8558, 
#    33.9785, -61.8171), cameraUpVector=(-0.682478, 0.730083, 0.0346623), 
#    cameraTarget=(-1.04391, -2.07928, 9.40767))
p = mdb.models['Beam'].parts['myBeam']
f = p.faces
p.Mirror(mirrorPlane=f[5], keepOriginal=OFF)
#session.viewports['Viewport: 1'].view.setValues(nearPlane=82.311, 
#    farPlane=156.478, width=79.1631, height=36.0925, cameraPosition=(86.103, 
#    80.4803, -18.6453), cameraUpVector=(-0.858047, 0.445569, 0.25539), 
#    cameraTarget=(-0.249799, -6.29825, 5.49083))
#session.viewports['Viewport: 1'].view.setValues(nearPlane=94.346, 
#    farPlane=168.431, width=90.7379, height=41.3697, cameraPosition=(83.2429, 
#    48.2362, 77.3385), cameraUpVector=(-0.407416, 0.699263, -0.587404), 
#    cameraTarget=(-0.12081, -4.84405, 1.16198))
#: Coordinates of vertex 2 :20.,15.,0.
#: Coordinates of vertex 1 :20.,15.,-20.
#session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
#    engineeringFeatures=ON)
#session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
#    referenceRepresentation=OFF)
mdb.models['Beam'].Material(name='Material-1')
mdb.models['Beam'].materials['Material-1'].Elastic(table=((1.6199424e+11, 0.28), ))
mdb.models['Beam'].HomogeneousSolidSection(name='Section-1', 
    material='Material-1', thickness=None)
p = mdb.models['Beam'].parts['myBeam']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
region = p.Set(cells=cells, name='Set-1')
p = mdb.models['Beam'].parts['myBeam']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
a = mdb.models['Beam'].rootAssembly
#session.viewports['Viewport: 1'].setValues(displayedObject=a)
#session.viewports['Viewport: 1'].assemblyDisplay.setValues(
#    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Beam'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Beam'].parts['myBeam']
a.Instance(name='myBeam-1', part=p, dependent=OFF)
#session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
#session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
#    meshTechnique=ON)
a = mdb.models['Beam'].rootAssembly
c1 = a.instances['myBeam-1'].cells
pickedRegions = c1.getSequenceFromMask(mask=('[#1 ]', ), )
a.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
elemType1 = mesh.ElemType(elemCode=C3D20R)
elemType2 = mesh.ElemType(elemCode=C3D15)
elemType3 = mesh.ElemType(elemCode=C3D10)
a = mdb.models['Beam'].rootAssembly
c1 = a.instances['myBeam-1'].cells
cells = c1.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(cells, )
a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
a = mdb.models['Beam'].rootAssembly
c1 = a.instances['myBeam-1'].cells
cells1 = c1.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(cells1, )
a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
a = mdb.models['Beam'].rootAssembly
e1 = a.instances['myBeam-1'].edges
pickedEdges2 = e1.getSequenceFromMask(mask=('[#2 ]', ), )

#myMinSize = 6
#myMaxSize = 20
myMinSize = 60
myMaxSize = 120

a.seedEdgeByBias(biasMethod=SINGLE, end2Edges=pickedEdges2, minSize=myMinSize, 
    maxSize=myMaxSize, constraint=FINER)
a = mdb.models['Beam'].rootAssembly
e1 = a.instances['myBeam-1'].edges
pickedEdges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
a.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, minSize=myMinSize, 
    maxSize=myMaxSize, constraint=FINER)
a = mdb.models['Beam'].rootAssembly
e1 = a.instances['myBeam-1'].edges
pickedEdges2 = e1.getSequenceFromMask(mask=('[#200 ]', ), )
a.seedEdgeByBias(biasMethod=SINGLE, end2Edges=pickedEdges2, minSize=myMinSize, 
    maxSize=myMaxSize, constraint=FINER)
a = mdb.models['Beam'].rootAssembly
e1 = a.instances['myBeam-1'].edges
pickedEdges1 = e1.getSequenceFromMask(mask=('[#200 ]', ), )
a.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, minSize=myMinSize, 
    maxSize=myMaxSize, constraint=FINER)
a = mdb.models['Beam'].rootAssembly
e1 = a.instances['myBeam-1'].edges
pickedEdges1 = e1.getSequenceFromMask(mask=('[#8 ]', ), )
a.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, minSize=myMinSize, 
    maxSize=myMaxSize, constraint=FINER)
#session.viewports['Viewport: 1'].view.setValues(nearPlane=87.4539, 
#    farPlane=161.334, width=83.8339, height=38.3476, cameraPosition=(52.8975, 
#    84.0614, 62.3706), cameraUpVector=(0.0478717, 0.292439, -0.955085), 
#    cameraTarget=(-0.915443, -1.98662, -10.2229))
#session.viewports['Viewport: 1'].view.setValues(nearPlane=90.7642, 
#    farPlane=160.759, width=87.0073, height=39.7991, cameraPosition=(83.2491, 
#    -22.6441, 79.3343), cameraUpVector=(-0.2236, 0.974342, -0.0257042), 
#    cameraTarget=(-1.00946, -1.6561, -10.2754))
#session.viewports['Viewport: 1'].view.setValues(nearPlane=93.8041, 
#    farPlane=154.896, width=89.9214, height=41.1321, cameraPosition=(-24.4754, 
#    60.5159, 96.0294), cameraUpVector=(0.600441, 0.548688, -0.581732), 
#    cameraTarget=(-1.85076, -1.00665, -10.145))
#session.viewports['Viewport: 1'].view.setValues(nearPlane=86.1752, 
#    farPlane=161.89, width=82.6083, height=37.7869, cameraPosition=(-89.0753, 
#    40.48, 68.7365), cameraUpVector=(0.715744, 0.691585, -0.0970568), 
#    cameraTarget=(-1.62776, -0.937486, -10.0508))
#session.viewports['Viewport: 1'].view.setValues(nearPlane=86.9039, 
#    farPlane=161.161, width=83.3069, height=38.1065, cameraUpVector=(0.556003, 
#    0.769129, -0.315121), cameraTarget=(-1.62776, -0.937485, -10.0508))
a = mdb.models['Beam'].rootAssembly
e1 = a.instances['myBeam-1'].edges
pickedEdges1 = e1.getSequenceFromMask(mask=('[#40 ]', ), )
a.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, minSize=myMinSize, 
    maxSize=myMaxSize, constraint=FINER)
a = mdb.models['Beam'].rootAssembly
partInstances =(a.instances['myBeam-1'], )
a.generateMesh(regions=partInstances)
#session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, 
#    adaptiveMeshConstraints=ON)
#session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
#    meshTechnique=OFF)
mdb.models['Beam'].StaticStep(name='Step-1', previous='Initial')
#session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
#session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
#    predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
a = mdb.models['Beam'].rootAssembly
s1 = a.instances['myBeam-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#20 ]', ), )
region = a.Surface(side1Faces=side1Faces1, name='Surf-1')
mdb.models['Beam'].Pressure(name='Load-1', createStepName='Step-1', 
    region=region, distributionType=UNIFORM, field='', magnitude=1e-06, 
    amplitude=UNSET)
#session.viewports['Viewport: 1'].view.setValues(nearPlane=85.8558, 
#    farPlane=161.228, width=82.5725, height=37.6469, cameraPosition=(-111.428, 
#    31.4839, -58.6459), cameraUpVector=(0.542434, 0.82198, 0.173535), 
#    cameraTarget=(-1.49313, -0.883298, -9.28354))
a = mdb.models['Beam'].rootAssembly
f1 = a.instances['myBeam-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#10 ]', ), )
region = a.Set(faces=faces1, name='Set-16')
mdb.models['Beam'].DisplacementBC(name='BC-1', createStepName='Step-1', 
    region=region, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
#session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
#    predefinedFields=OFF, connectors=OFF)
