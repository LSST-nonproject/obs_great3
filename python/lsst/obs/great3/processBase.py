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
from lsst.meas.base import SingleFrameMeasurementTask

from .buildPsf import *

class ProcessBaseConfig(lsst.pex.config.Config):

    psf = lsst.pex.config.ConfigurableField(
        target = BuildControlPsfTask,
        doc = "PSF determination"
    )
    measurement = lsst.pex.config.ConfigurableField(
        target = SingleFrameMeasurementTask,
        doc = "Measurement of source properties"
    )
    varianceBorderWidth = lsst.pex.config.Field(
        dtype=int,
        default=1,
        doc=("Use a border of this many pixels around each postage stamp (combined over the full image)"
             " to compute the variance")
    )
    galaxyStampSize = lsst.pex.config.Field(
        dtype=int,
        default=48,
        optional=False,
        doc="Cutout size for galaxies"
    )

    def setDefaults(self):
        lsst.pex.config.Config.setDefaults(self)
        self.measurement.slots.centroid = "base_SdssCentroid"
        self.measurement.slots.instFlux = None
        self.measurement.slots.modelFlux = None
        self.measurement.slots.apFlux = "base_GaussianFlux"
        self.measurement.algorithms.names = ["base_SdssCentroid", "base_SdssShape", "base_GaussianFlux", "base_PsfFlux"]

class ProcessBaseTask(lsst.pipe.base.CmdLineTask):

    ConfigClass = ProcessBaseConfig

    def __init__(self, **kwds):
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        self.algMetadata = lsst.daf.base.PropertyList()
        self.makeSubtask("psf")
        self.makeSubtask("measurement", schema=self.schema, algMetadata=self.algMetadata)

    def computeVariance(self, image):
        array = image.getArray()
        nGals = array.shape[0] / self.config.galaxyStampSize
        assert nGals * self.config.galaxyStampSize == array.shape[0]
        n = array.shape[0] / nGals
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
