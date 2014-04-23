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

class ReductionTaskRunner(lsst.pipe.base.TaskRunner):
    """A TaskRunner that passes all DataRefs to the run() method at once.
    """

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(parsedCmd.butler, parsedCmd.id.refList,)]

class MakeShearCatalogConfig(lsst.pex.config.Config):
    useWeights = lsst.pex.config.Field(
        dtype=bool,
        doc="Whether to use per-object SNR weights in shear estimation",
        default=True,
        )
    distortionField = lsst.pex.conofig.Field(
        dtype=str,
        doc=("Prefix for field containing distortion values in .e1, .e2, and error in .err"),
        default="shape.hsm.regauss",
        )
    useImprovedSigmaFactor = lsst.pex.config.Field(
        dtype=bool,
        doc="Whether to use the improved 2/(1+2*sigma_g^2) factor from BJ02 in converting sigma_g to sigma_e"
        default=True,
        )
    useImprovedResponsivity = lsst.pex.config.Field(
        dtype=bool,
        doc="Whether to use the improved responsivity factor from BJ02 (eqn 5-33) in distortion to shear"
        default=False,
        )
    shearCalibrationFactor = lsst.pex.config.Field(
        dtype=float,
        doc="Empirical shear calibration factor; default value is appropriate for Regaussianization",
        default=0.98
        )
    nRmsIterations = lsst.pex.config.Field(
        dtype=int,
        doc=("Number of iterations between estimating weights from RMS and estimating RMS from weights "
             "(i.e. 0 to use unweighted RMS estimate)"),
        default=1,
        )

class MakeShearCatalogTask(lsst.pipe.base.CmdLineTask):

    RunnerClass = ReductionRunner
    ConfigClass = MakeShearCatalogConfig

    _DefaultName = "makeShearCatalog"

    def __init__(self, **kwds):
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        self.schema = lsst.afw.table.Schema()
        self.idKey = self.schema.addField("id", type=long, doc="unique ID")
        self.e1Key = self.schema.addField("e1", type=float, doc="axis-aligned distortion component")
        self.e2Key = self.schema.addField("e2", type=float, doc="off-axis distortion component")
        self.errKey = self.schema.addField("err", type=float, doc="raw error ('sigma_g')")
        self.eSigmaKey = self.schema.addField("eSigma", type=float,
                                              doc="corrected error on distortion ('sigma_e')")
        self.g1Key = self.schema.addField("g1", type=float, doc="axis-aligned shear component")
        self.g2Key = self.schema.addField("g2", type=float, doc="off-axis shear component")
        if self.config.useWeights:
            self.weightKey = self.schema.addField("weight", type=float,
                                                  doc="per-object weight to apply when combining shears")

    def run(self, butler, dataRefList):
        # Combine the catalogs from the different subfields, copying only the relevant values
        mapper = None
        outCat = lsst.afw.table.BaseCatalog(self.schema)
        outCat.reserve(10000*len(dataRefList))
        for dataRef in dataRefList:
            self.log.info("Reading catalog for subfield %d" % dataRef.dataId["subfield"])
            inCat = dataRef.get("src", immediate=True)
            if mapper is None:
                flagKey = inCat.schema.find(self.config.distortionField + ".flags").key
                mapper = lsst.afw.table.SchemaMapper(inCat.schema, outCat.schema)
                mapper.addMapping(inCat.schema.find("id").key, "id", True)
                mapper.addMapping(inCat.schema.find(self.config.distortionField + ".e1").key, "e1", True))
                mapper.addMapping(inCat.schema.find(self.config.distortionField + ".e2").key, "e2", True))
                mapper.addMapping(inCat.schema.find(self.config.distortionField + ".err").key, "err", True))
                assert mapper.getOutputSchema() == self.schema
            outCat.extend(inCat[numpy.logical_not(inCat.get(flagKey))], mapper=mapper)

        # Correct uncertainties to be applicable to distortion, not shear
        self.log.info("Computing per-object distortion uncertainty")
        if self.config.useImproveSigmaFactor:
            factor = 2.0 / (1.0 + 2.0*outCat[self.errKey]**2)
        else:
            factor = 2.0
        outCat[self.eSigmaKey][:] = factor * outCat[self.errKey]

        # Iteratively compute RMS weights and ellipticity (note that default is just one iteration)
        def computeRMS(e1, e2, eSigma, w):
            return numpy.sum(0.5*(e1**2 + e2**2 - 2*eSigma)*w)
        rms = computeRMS(outCat[self.e1Key], outCat[self.e2Key], outCat[self.eSigmaKey], 1.0 / len(outCat))
        self.log.info("Unweighted distortion RMS=%s" % rms)
        if self.config.useWeights:
            outCat[self.weightKey][:] = 1.0 / (rms**2 + outCat[self.eSigmaKey]**2)
            outCat[self.weightKey][:] /= numpy.sum(outCat[self.weightKey])
            for i in range(self.config.nRmsIterations):
                rms = computeRMS(outCat[self.e1Key], outCat[self.e2Key], outCat[self.eSigmaKey],
                                 outCat[self.weightKey])
                self.log.info("Iteration %d, weighted distortion RMS=%s" % (i, rms))

        # Compute shears, including responsivity
        if self.config.useImprovedResponsivity:
            raise NotImplementedError("Improved responsivity is not yet implemented")
        else:
            responsivity = 2.0 * (1.0 - rms**2)
            if responsivity > 2.0 or responsivity < 1.0:
                raise ValueError("Invalid responsivity: %s" % responsivity)
        self.log.info("Computing shears with calibration factor %s and responsivity %s"
                      % (self.config.shearCalibrationFactor, responsivity))
        f = self.config.shearCalibrationFactor / responsivity
        outCat[self.g1Key][:] = f*outCat[self.e1Key]
        outCat[self.g2Key][:] = f*outCat[self.e2Key]

        # Write outputs
        butler.put(outCat, "shear")

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
