setImageType('FLUORESCENCE');
setPixelSizeMicrons(0.138000, 0.138000)
createSelectAllObject(true)
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage":"Channel 1","requestedPixelSizeMicrons":0.25,"backgroundRadiusMicrons":8.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.7,"minAreaMicrons":10.0,"maxAreaMicrons":400.0,"threshold":100.0,"watershedPostProcess":true,"cellExpansionMicrons":0.0,"includeNuclei":false,"smoothBoundaries":true,"makeMeasurements":true}')

import qupath.imagej.gui.IJExtension

double downsample = 1.0

def server = getCurrentServer()
def selectedObject = getSelectedObject()
def request = RegionRequest.createInstance(server.getPath(), downsample, selectedObject.getROI())
boolean setROI = true

IJExtension.getImageJInstance()

def imp = IJExtension.extractROIWithOverlay(
    getCurrentServer(),
    selectedObject,
    getCurrentHierarchy(),
    request,
    setROI,
    getCurrentViewer().getOverlayOptions()
    ).getImage()
 
 imp.show()

