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

from .buildPsf import *

class ProcessSingleEpochConfig(lsst.pex.config.Config):

    psf = lsst.pex.config.ConfigurableField(
        target = BuildControlPsfTask,
        doc = "PSF determination"
    )
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

    def setDefaults(self):
        lsst.pex.config.Config.setDefaults(self)
        self.measurement.slots.centroid = "centroid.sdss"
        self.measurement.slots.instFlux = None
        self.measurement.slots.modelFlux = None
        self.measurement.slots.apFlux = "flux.sinc"
        self.measurement.algorithms = ["shape.sdss", "flux.sinc", "flux.psf"]

class ProcessSingleEpochTask(lsst.pipe.base.CmdLineTask):

    ConfigClass = ProcessSingleEpochConfig

    _DefaultName = "ProcessSingleEpoch"

    def __init__(self, **kwds):
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        self.algMetadata = lsst.daf.base.PropertyList()
        self.makeSubtask("psf")
        self.makeSubtask("measurement", schema=self.schema, algMetadata=self.algMetadata)

    def computeVariance(self, image):
        array = image.getArray()
        n = array.shape[0] / 100
        assert n * 100 == array.shape[0]
        mask = numpy.zeros(array.shape, dtype=bool)
        for i in range(self.config.varianceBorderWidth):
            mask[i::n,:] = True
            mask[:,i::n] = True
            mask[n-i::n,:] = True
            mask[:,n-i::n] = True
        borderPixels = array[mask]
        return numpy.std(borderPixels, dtype=numpy.float64)**2

    def buildExposure(self, dataRef):
        image = dataRef.get("image", immediate=True)
        exposure = lsst.afw.image.ExposureF(image.getBBox(lsst.afw.image.PARENT))
        exposure.getMaskedImage().getImage().getArray()[:,:] = image.getArray()
        exposure.getMaskedImage().getVariance().set(self.computeVariance(image))
        exposure.setPsf(self.psf.run(dataRef))
        return exposure

    def buildSourceCatalog(self, imageBBox, dataRef):
        """Build an empty source catalog, using the provided sim catalog's position to generate
        square Footprints and its ID to set that of the source catalog.
        """
        sourceCat = lsst.afw.table.SourceCatalog(self.schema)
        sourceCat.getTable().setMetadata(self.algMetadata)
        simCat = dataRef.get("galaxy_catalog", immediate=True)
        xKey = simCat.schema.find('x').key
        yKey = simCat.schema.find('y').key
        idKey = simCat.schema.find('ID').key
        n = imageBBox.getWidth() / 100
        assert n * 100 == imageBBox.getWidth()
        dims = lsst.afw.geom.Extent2I(n, n)
        offset = lsst.afw.geom.Extent2I(simCat[0][xKey], simCat[0][yKey])
        for simRecord in simCat:
            sourceRecord = sourceCat.addNew()
            sourceRecord.setId(simRecord.get(idKey))
            position = lsst.afw.geom.Point2I(simRecord.get(xKey), simRecord.get(yKey))
            bbox = lsst.afw.geom.Box2I(position - offset, dims)
            footprint = lsst.afw.detection.Footprint(bbox, imageBBox)
            footprint.getPeaks().push_back(lsst.afw.detection.Peak(position.getX(), position.getY()))
            sourceRecord.setFootprint(footprint)
        return sourceCat

    def run(self, dataRef):
        exposure = self.buildExposure(dataRef)
        sourceCat = self.buildSourceCatalog(exposure.getBBox(lsst.afw.image.PARENT), dataRef)
        self.measurement.run(exposure, sourceCat)
        dataRef.put(sourceCat, "src")

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="image", level="image",
                               help="data ID, e.g. --id subfield=0")
        return parser

    def writeConfig(self, butler, clobber=False):
        pass

    def writeSchemas(self, butler, clobber=False):
        pass

    def writeMetadata(self, dataRef):
        pass

    def writeEupsVersions(self, butler, clobber=False):
        pass
