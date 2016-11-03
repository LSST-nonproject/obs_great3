#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 20082014 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import numpy
import lsst.daf.base
import lsst.pex.config
import lsst.pipe.base
import lsst.afw.table
import lsst.afw.image
from lsst.meas.algorithms import  SourceDetectionTask, SourceMeasurementTask
from lsst.meas.deblender import SourceDeblendTask
try:
    import scipy.spatial
    spatialAvailable = True
except ImportError:
    spatialAvailable = False



from .processBase import *

class ProcessSingleEpochConfig(ProcessBaseConfig):
    doDetection = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc=("run source detection on the image")
    )
    detection = lsst.pex.config.ConfigurableField(
        target = SourceDetectionTask,
        doc = "detection algorithm",
    )
    deblend = lsst.pex.config.ConfigurableField(
        target = SourceDeblendTask,
        doc = "deblend algorithm",
    )
    addGridPoints = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc=("add footprints at the grid points if they are not in the detection list")
    )
    addGridDist = lsst.pex.config.Field(
        dtype=float,
        default=6.,
        doc=("add objects grid points that do not have a source within this distance")
    )
    correlatedNoiseFile = lsst.pex.config.Field(
        doc = "File that contains list of k values and power spectrum to use with correlated noise.  Still" \
        " needs to be scaled by noise variance.  For BFD only",
        dtype = str,
        default = '',
    )
    doCorrelatedNoise = lsst.pex.config.Field(
        doc = "Run with correlated noise",
        dtype = bool,
        default = False,
    )
    recomputeVariance = lsst.pex.config.Field(
        doc = "recompute Variance after deblending",
        dtype = bool,
        default = False,
    )

    def setDefaults(self):
        ProcessBaseConfig.setDefaults(self)
        self.detection.thresholdValue = 5
        self.detection.reEstimateBackground = False


class ProcessSingleEpochTask(ProcessBaseTask):

    ConfigClass = ProcessSingleEpochConfig
    _DefaultName = "processSingleEpoch"

    def __init__(self, **kwargs):
        super(ProcessSingleEpochTask, self).__init__(**kwargs)

        if self.config.doDetection:

            self.makeSubtask("detection", schema=self.schema)
            self.makeSubtask("deblend", schema=self.schema)

            self.indexKey = self.schema.addField("index", type=int, doc="grid index of galaxy")
            self.gridDistKey = self.schema.addField("gridDist", type=float, doc="distance to nearest grid point")

        self.gridPositions = []
        self.gridBoxes = []

    def buildSourceCatalog(self, imageBBox, dataRef):
        """Build an empty source catalog, using the provided sim catalog's position to generate
        square Footprints and its ID to set that of the source catalog.
        """
        sourceCat = lsst.afw.table.SourceCatalog(self.schema)
        sourceCat.getTable().setMetadata(self.algMetadata)
        simCat = dataRef.get(self.config.dataType + "epoch_catalog", immediate=True)
        xKey = simCat.schema.find('x').key
        yKey = simCat.schema.find('y').key
        idKey = simCat.schema.find('num').key
        n = imageBBox.getWidth() / self.config.numPerRow
        assert n * self.config.numPerRow == imageBBox.getWidth()
        dims = lsst.afw.geom.Extent2I(n, n)
        offset = lsst.afw.geom.Extent2I(int(numpy.min(simCat[xKey])), int(numpy.min(simCat[yKey])))

        if self.config.maxObjects is not None:
            simCat = simCat[:self.config.maxObjects]

        for simRecord in simCat:
            sourceRecord = sourceCat.addNew()
            sourceRecord.setId(simRecord.get(idKey))
            position = lsst.afw.geom.Point2I(int(simRecord.get(xKey)), int(simRecord.get(yKey)))
            bbox = lsst.afw.geom.Box2I(position - offset, dims)
            footprint = lsst.afw.detection.Footprint(bbox, imageBBox)
            # add dummy value of 100 for peak value
            footprint.addPeak(position.getX(), position.getY(), 100)
            self.gridPositions.append((position.getX(),position.getY()))
            self.gridBoxes.append(bbox)
            sourceRecord.setFootprint(footprint)

        if self.config.doCorrelatedNoise and bfdAvailable:
            data = numpy.genfromtxt(self.config.correlatedNoiseFile)
            sourceRecord.getTable().getMetadata().set('kData',data[:,0])
            sourceRecord.getTable().getMetadata().set('psData',data[:,1])

        return sourceCat

    def run(self, dataRef):
        exposure = self.buildExposure(dataRef)
        sourceCat = self.buildSourceCatalog(exposure.getBBox(lsst.afw.image.PARENT), dataRef)

        if self.config.doDetection:
            table = lsst.afw.table.SourceTable.make(self.schema)
            table.setMetadata(self.algMetadata)

            detections = self.detection.makeSourceCatalog(table, exposure)
            sourceCat = detections.sources
            if self.config.maxObjects is not None:
                sourceCat = sourceCat[:self.config.maxObjects]

            if self.config.addGridPoints:
                if spatialAvailable:
                    lastId = sourceCat[-1].getId() + 1

                    peakPosList = [(peak.getIx(), peak.getIy()) for src in sourceCat
                                   for peak in src.getFootprint().getPeaks()]

                    gridPos = numpy.array(self.gridPositions)
                    peakPos = numpy.array(peakPosList)
                    tree = scipy.spatial.cKDTree(peakPos)

                    for (x,y),box in zip(self.gridPositions, self.gridBoxes):
                        dist, index = tree.query([x,y])
                        gridDetection = False
                        if dist < self.config.addGridDist:
                            gridDetection = True

                        if gridDetection is False:
                            sourceRecord = sourceCat.addNew()
                            sourceRecord.setId(lastId+1)
                            lastId += 1
                            footprint = lsst.afw.detection.Footprint(box, exposure.getBBox())
                            footprint.addPeak(x, y, 100)
                            sourceRecord.setFootprint(footprint)
                else:
                    self.log.warn('Cannot add grid points requires scipy.spatial')

            if self.config.recomputeVariance:
                mask = exposure.getMaskedImage().getMask()
                mask &= (mask.getPlaneBitMask("DETECTED"))
                iso_mask = (mask.getArray() != (mask.getPlaneBitMask("DETECTED")))
                variance = numpy.var(exposure.getMaskedImage().getImage().getArray()[iso_mask], dtype=numpy.float64)
                self.log.info('Computed variance after detection: %f' % variance)
                exposure.getMaskedImage().getVariance().set(variance)
                self.algMetadata.set('noise_variance',variance)

            self.deblend.run(exposure, sourceCat, exposure.getPsf())

        self.measurement.run(exposure, sourceCat)

        if self.config.doDetection:

            # Remove parents from catalog, only keeping the children
            mask = numpy.array([a.get('deblend.nchild') == 0 for a in sourceCat])
            sourceCat = sourceCat.subset(mask)

            if spatialAvailable:
                gridPos = numpy.array(self.gridPositions)
                tree = scipy.spatial.cKDTree(gridPos)

                srcPosList = numpy.array([(src.getX(), src.getY()) for src in sourceCat])
                minDist, index = tree.query(srcPosList)

                # assume a square image
                stampSize = exposure.getMaskedImage().getImage().getWidth() / self.config.numPerRow
                for dist,src in zip(minDist, sourceCat):

                    xIndex = src.getX() // stampSize
                    yIndex = src.getY() // stampSize
                    index = xIndex*self.config.numPerRow + yIndex
                    src.set(self.gridDistKey, dist)
                    src.set(self.indexKey, int(index))
            else:
                self.log.warn("Not computing distances or indexes because can't find scipy.spatial")


        dataRef.put(sourceCat, self.config.dataType + "src")

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="image", level="image",
                               help="data ID, e.g. --id subfield=0")
        return parser

    def writeConfig(self, butler, clobber=False, doBackup=False):
        pass

    def writeSchemas(self, butler, clobber=False, doBackup=False):
        pass

    def writeMetadata(self, dataRef):
        pass

    def writeEupsVersions(self, butler, clobber=False, doBackup=False):
        pass
