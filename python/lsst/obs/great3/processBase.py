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
import lsst.afw.coord
import lsst.afw.geom
import lsst.meas.extensions.shapeHSM
from lsst.meas.algorithms import SourceMeasurementTask
import lsst.meas.multifit
import lsst.meas.extensions.multiShapelet
import lsst.meas.extensions.photometryKron
import lsst.obs.great3
try:
    import lsst.meas.extensions.bfd
    bfdAvailable = True
except ImportError:
    bfdAvailable = False
from .buildPsf import *

class ProcessBaseConfig(lsst.pex.config.Config):

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
    dataType = lsst.pex.config.Field(
        dtype=str,
        default='',
        doc=("default data type")
    )
    maxObjects = lsst.pex.config.Field(
        dtype=int,
        default=None,
        optional=True,
        doc=("If this is not None, then only process this many objects.")
    )

    def setDefaults(self):
        lsst.pex.config.Config.setDefaults(self)
        self.measurement.slots.centroid = "centroid.sdss"
        self.measurement.slots.instFlux = None
        self.measurement.slots.modelFlux = None
        self.measurement.slots.calibFlux = None
        self.measurement.slots.apFlux = "flux.sinc"
        self.measurement.algorithms = ["shape.sdss", "flux.sinc", "flux.psf", "shape.hsm.regauss", "flux.kron",
                                       "cmodel", "classification.extendedness"]
        self.measurement.algorithms.names |= lsst.meas.extensions.multiShapelet.algorithms
        self.measurement.algorithms['classification.extendedness'].fluxRatio = 0.985
        if bfdAvailable:
            self.measurement.algorithms.names |= ["bfd.kmoment"]
            self.measurement.algorithms['bfd.kmoment'].sigma = 2
            self.measurement.algorithms['bfd.kmoment'].useRecVariance=False
            self.measurement.algorithms['bfd.kmoment'].useTableVariance=True
            self.measurement.algorithms['bfd.kmoment'].shift = True
            self.measurement.algorithms['bfd.kmoment'].wIndex = 4


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

        image = dataRef.get(self.config.dataType + "image", immediate=True)
        exposure = lsst.afw.image.ExposureF(image.getBBox(lsst.afw.image.PARENT))
        exposure.getMaskedImage().getImage().getArray()[:,:] = image.getArray()

        variance = self.computeVariance(image)
        exposure.getMaskedImage().getVariance().set(variance)
        self.algMetadata.set('noise_variance',variance)

        exposure.setPsf(self.psf.run(dataRef, self.config.dataType))

        pixelScale = dataRef.getButler().mapper.PIXEL_SCALE
        wcs = lsst.afw.image.makeWcs(lsst.afw.coord.Coord(0.*lsst.afw.geom.degrees, 0.*lsst.afw.geom.degrees),
                                     lsst.afw.geom.Point2D(0.0,0.0),
                                     pixelScale.asDegrees(), 0.0, 0.0,
                                     pixelScale.asDegrees())
        exposure.setWcs(wcs)

        calib = lsst.afw.image.Calib()
        calib.setFluxMag0(1e12)
        exposure.setCalib(calib)

        return exposure
