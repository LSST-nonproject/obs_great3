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
from lsst.meas.algorithms import SourceMeasurementTask

class OutputOnlyDataIdContainer(lsst.pipe.base.DataIdContainer):
    """A replacement for the basic DataIdContainer that doesn't demand that its
    DataRefs pass a datasetExists test.
    """

    def makeDataRefList(self, namespace):
        """Compute refList based on idList

        Not called if add_id_argument called with doMakeDataRef=False

        @param namespace: results of parsing command-line (with 'butler' and 'log' elements)
        """
        if self.datasetType is None:
            raise RuntimeError("Must call setDatasetType first")
        butler = namespace.butler
        for dataId in self.idList:
            self.refList += list(butler.subset(datasetType=self.datasetType, level=self.level, dataId=dataId))



class TileStarsConfig(lsst.pex.config.Config):
    nSubtileX = lsst.pex.config.Field(dtype=int, default=10,
                                      doc="Number of partitions of each tile in x direction")
    nSubtileY = lsst.pex.config.Field(dtype=int, default=10,
                                      doc="Number of partitions of each tile in y direction")
    measurement = lsst.pex.config.ConfigurableField(
        target = SourceMeasurementTask,
        doc = "Measurement of source properties"
    )
    varianceBorderWidth = lsst.pex.config.Field(
        dtype=int,
        default=1,
        doc=("Use a border of this many pixels around each postage stamp (combined over the full image)"
             " to compute the variance")
    )
    doSkipCompleted = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Skip data IDs that have already been processed"
    )

    def setDefaults(self):
        lsst.pex.config.Config.setDefaults(self)
        self.measurement.centroider = "centroid.gaussian"
        self.measurement.slots.centroid = "centroid.gaussian"
        self.measurement.slots.instFlux = None
        self.measurement.slots.modelFlux = None
        self.measurement.slots.apFlux = "flux.sinc"
        self.measurement.algorithms = ["shape.sdss", "flux.sinc", "flux.psf"]

class TileStarsTask(lsst.pipe.base.CmdLineTask):

    ConfigClass = TileStarsConfig

    _DefaultName = "tileStars"

    def __init__(self, **kwds):
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        self.algMetadata = lsst.daf.base.PropertyList()
        self.makeSubtask("measurement", schema=self.schema, algMetadata=self.algMetadata)

    def readIndex(self, dataRef):
        TILE_SIZE = dataRef.butlerSubset.butler.mapper.TILE_SIZE
        PIXEL_SCALE = dataRef.butlerSubset.butler.mapper.PIXEL_SCALE
        tilePixelWidth = int(numpy.rint(TILE_SIZE.asDegrees() / PIXEL_SCALE.asDegrees()))
        tilePixelHeight = tilePixelWidth
        subtilePixelWidth = tilePixelWidth // self.config.nSubtileX
        subtilePixelHeight = tilePixelHeight // self.config.nSubtileY
        if not tilePixelWidth % self.config.nSubtileX == 0:
            raise ValueError("Tile width in pixels %d does not divide evenly into %d segments"
                             % (tilePixelWidth, self.config.nSubtileX))
        if not tilePixelHeight % self.config.nSubtileY == 0:
            raise ValueError("Tile height in pixels %d does not divide evenly into %d segments"
                             % (tilePixelHeight, self.config.nSubtileY))
        fullBBox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(dataRef.dataId['sx'] * subtilePixelWidth,
                                                             dataRef.dataId['sy'] * subtilePixelHeight),
                                       lsst.afw.geom.Extent2I(subtilePixelWidth, subtilePixelHeight))
        index = dataRef.get("star_index", immediate=True)
        conditions = [
            index['tile.index.x'] == dataRef.dataId['tx'],
            index['tile.index.y'] == dataRef.dataId['ty'],
            index['tile.position.x'] > fullBBox.getMinX(),
            index['tile.position.y'] > fullBBox.getMinY(),
            index['tile.position.x'] < fullBBox.getMaxX(),
            index['tile.position.y'] < fullBBox.getMaxY()
            ]
        index[numpy.logical_and.reduce(conditions)]
        return index, fullBBox

    def build(self, dataRef, index, fullBBox):
        """Build the exposure and a mostly-empty SourceCatalog containing the objects we put into it.
        """
        catalog = lsst.afw.table.SourceCatalog(self.schema)
        PSTAMP_SIZE = dataRef.butlerSubset.butler.mapper.PSTAMP_SIZE
        exposure = lsst.afw.image.ExposureF(fullBBox)
        mi = exposure.getMaskedImage()
        mask = mi.getMask()
        mask.addMaskPlane("VAREST") # used to set a border region we can use to compute the variance
        detectedBitMask = mask.getPlaneBitMask("DETECTED")
        varEstBitMask = mask.getPlaneBitMask("VAREST")

        xInKey = index.schema.find('subfield.position.x').key
        yInKey = index.schema.find('subfield.position.y').key
        xOutKey = index.schema.find('tile.position.x').key
        yOutKey = index.schema.find('tile.position.y').key
        subfieldKey = index.schema.find('subfield.index').key
        pStampDim = lsst.afw.geom.Extent2I(PSTAMP_SIZE, PSTAMP_SIZE)
        pStampOffset = lsst.afw.geom.Extent2I(index[0].get(xInKey) % PSTAMP_SIZE,
                                              index[0].get(yInKey) % PSTAMP_SIZE)
        varEstStencil = lsst.afw.image.MaskU(pStampDim)
        varEstStencil.getArray()[:self.config.varianceBorderWidth,:] = varEstBitMask
        varEstStencil.getArray()[-self.config.varianceBorderWidth:,:] = varEstBitMask
        varEstStencil.getArray()[:,:self.config.varianceBorderWidth] = varEstBitMask
        varEstStencil.getArray()[:,-self.config.varianceBorderWidth:] = varEstBitMask

        detectedStencil = lsst.afw.image.MaskU(pStampDim)
        detectedStencil.set(detectedBitMask)

        for record in index:
            outCenter = lsst.afw.geom.Point2I(record.get(xOutKey), record.get(yOutKey))
            outBBox = lsst.afw.geom.Box2I(outCenter - pStampOffset, pStampDim)
            if not fullBBox.contains(outBBox):  # skip partial overlaps
                self.log.logdebug("Skipping star at %d,%d on edge" % (outCenter.getX(), outCenter.getY()))
                continue
            outSub = mi.Factory(mi, outBBox, lsst.afw.image.PARENT)
            outSubImage = outSub.getImage()
            outSubMask = outSub.getMask()
            if (outSubMask.getArray() & detectedBitMask).any():
                self.log.logdebug("Skipping star at %d,%d due to overlap with existing star"
                                  % (outCenter.getX(), outCenter.getY()))
                continue
            outSubMask |= varEstStencil
            outSubMask |= detectedStencil

            inCenter = lsst.afw.geom.Point2I(record.get(xInKey), record.get(yInKey))
            inBBox = lsst.afw.geom.Box2I(inCenter - pStampOffset, pStampDim)

            inSubImage = dataRef.get("starfield_image_sub", subfield=record.get(subfieldKey), bbox=inBBox,
                                     immediate=True)
            outSubImage += inSubImage

            source = catalog.addNew()
            footprint = lsst.afw.detection.Footprint(outBBox)
            footprint.getPeaks().push_back(lsst.afw.detection.Peak(outCenter.getX(), outCenter.getY()))
            source.setFootprint(footprint)

        varEstPixels = (mi.getMask().getArray() & varEstBitMask).astype(bool)
        emptyPixels = numpy.logical_not((mi.getMask().getArray() & detectedBitMask).astype(bool))
        stddev = numpy.std(mi.getImage().getArray()[varEstPixels], dtype=numpy.float64)
        mi.getImage().getArray()[emptyPixels] = numpy.random.randn(emptyPixels.sum()) * stddev
        mi.getVariance().set(stddev**2)
        mask.removeAndClearMaskPlane("VAREST")

        return exposure, catalog

    def run(self, dataRef):
        if self.config.doSkipCompleted and dataRef.datasetExists("subtile_star_catalog"):
            self.log.info("Skipping %s; already completed" % dataRef.dataId)

        self.log.info("Tiling stars for %s" % dataRef.dataId)

        seed = hash(tuple(dataRef.dataId.values()))
        if seed < 0: seed *= -1
        numpy.random.seed(seed)

        index, fullBBox = self.readIndex(dataRef)

        exposure, catalog = self.build(dataRef, index, fullBBox)

        self.measurement.run(exposure, catalog)

        dataRef.put(exposure, "subtile_star_image")
        dataRef.put(catalog, "subtile_star_catalog")

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="subtile_star_image",
                               ContainerClass=OutputOnlyDataIdContainer,
                               help="data ID, e.g. --id field=0 epoch=0 tx=0..4 ty=0..4 sx=0..10 sy=0..10")
        return parser

    def writeConfig(self, butler, clobber=False):
        pass

    def writeSchemas(self, butler, clobber=False):
        pass

    def writeMetadata(self, dataRef):
        pass

    def writeEupsVersions(self, butler, clobber=False):
        pass
