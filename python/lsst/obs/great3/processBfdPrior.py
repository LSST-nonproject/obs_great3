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
import lsst.meas.extensions.bfd

from .processSingleEpoch import *

class ProcessBfdPriorConfig(ProcessSingleEpochConfig):
    snMin = lsst.pex.config.Field(
        dtype=float,
        default=5,
        optional=True,
        doc="Minimun flux S/N"
    )
    snMax = lsst.pex.config.Field(
        dtype=float,
        default=40,
        optional=True,
        doc="Maximum flux S/N"
    )
    noiseFactor = lsst.pex.config.Field(
        dtype=float,
        default=1,
        optional=True,
        doc="Noise boost factor for kernel smoothing"
    )
    priorSigmaCutoff = lsst.pex.config.Field(
        dtype=float,
        default=5.5,
        optional=True,
        doc="Maximum sigma range when sampling for prior"
    )
    priorSigmaStep = lsst.pex.config.Field(
        dtype=float,
        default=1.,
        optional=True,
        doc="Step size when sampling for prior"
    )
    priorSigmaBuffer = lsst.pex.config.Field(
        dtype=float,
        default=1.,
        optional=True,
        doc="Buffer width of KdTreePrior (in sigma)"
    )
    nSample = lsst.pex.config.Field(
        dtype=int,
        default=30000,
        optional=True,
        doc="Number of templates sampled per target"
    )
    maxXY = lsst.pex.config.Field(
        dtype=float,
        default=4.,
        optional=True,
        doc="Maximum translational displacement in sigma of the nominal covariance matrix"
    )
    sigma = lsst.pex.config.Field(
        dtype=float,
        default=2,
        optional=True,
        doc="Sigma used in k-space weight function"
    )
    wIndex = lsst.pex.config.Field(
        dtype=int,
        default=4,
        optional=True,
        doc="index used in k-space weight function"
    )
    centroidName = lsst.pex.config.Field(
        dtype=str,
        default='centroid.sdss',
        optional=True,
        doc="name of centroid to use from the catalog"
    )
    selectionOnly = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="only do selection in the prior"
    )
    covFile = lsst.pex.config.Field(
        dtype=str,
        default='./test.fits',
        optional=True,
        doc="file that contains the covariance matrices"
    )
    covFromCat = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="get covariance matrix from the first entry in a catalog"
    )
    isoCov = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        optional=True,
        doc="isotropize covariance matrix"
    )
    sample = lsst.pex.config.Field(
        dtype=float,
        default=0.2,
        optional=True,
        doc="Only use this fraction of the galaxies"
    )
    sampleDict = lsst.pex.config.DictField(
        keytype=str, itemtype=float,
        default={},
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for each covariance label"
    )
    useLabels = lsst.pex.config.ListField(
        dtype=str,
        default=[],
        optional=True,
        doc="List of labels for which to build the prior. "
    )
    invariantCovariance = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="will all the galaxies have the same covariance"
    )
    reCentroid = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="recentroid galaxis"
    )
    reCentroidPsf = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="recentroid galaxis"
    )
    scaleFactor = lsst.pex.config.Field(
        dtype=float,
        default=1./10,
        optional=True,
        doc="rescale images by this factor"
    )

    def setDefaults(self):
        self.dataType = 'deep_'

class ProcessBfdPriorTask(ProcessSingleEpochTask):

    ConfigClass = ProcessBfdPriorConfig
    _DefaultName = "processBfdPrior"


    def __init__(self,  **kwargs):
        """Initialize the measurement task, including the modelfits catalog schema,
        the model, prior, and calib objects, and the fitter subtask.
        """

        ProcessSingleEpochTask.__init__(self, **kwargs)
        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()

    def getCovariances(self):
        if not self.config.covFromCat:
            cat = lsst.afw.table.BaseCatalog.readFits(self.config.covFile)
            covList = []
            for rec in cat:
                cov = rec.get('isoCov')
                label = rec.get('label')
                minVariance = rec.get('min')
                maxVariance = rec.get('max')
                cov = numpy.array(cov.reshape(6,6), dtype=numpy.float32)

                covList.append((cov,label,minVariance,maxVariance))

            return covList
        else:
             cat = lsst.afw.table.BaseCatalog.readFits(self.config.covFile)
             mask = cat['bfd.flags'] == 0
             cov = numpy.array(cat[mask][0].get('bfd.momentsCov').reshape(6,6), dtype=numpy.float32)

             # make it isotropic
             if self.config.isoCov:
                 varE = 0.5*(cov[4,4] + cov[5,5])
                 varC = 0.5*(cov[1,1] + cov[2,2])
                 cov[1:3,:]=0
                 cov[:,1:3]=0
                 cov[4:6,:]=0
                 cov[:,4:6]=0
                 cov[1,1]=varC
                 cov[2,2]=varC
                 cov[4,4]=varE
                 cov[5,5]=varE

             return [(cov,'b0',-1,-1)]

    def run(self, dataRef):

        exposure = self.buildExposure(dataRef)
        sourceCat = dataRef.get('deep_src')
        covList = self.getCovariances()

        exposure.getMaskedImage().getImage().getArray()[:] *= self.config.scaleFactor
        exposure.getMaskedImage().getVariance().getArray()[:] *= self.config.scaleFactor**2

        for iter,(cov, label, minVar, maxVar) in enumerate(covList):

            if (label not in self.config.useLabels) and len(self.config.useLabels) > 0:
                continue

            if self.config.sampleDict is None:
                sample = self.config.sample
            else:
                if label not in self.config.sampleDict.keys():
                    self.log.info('Cannot find label %s, in sampleDict will use the default value %f' %
                                  (label,self.config.sample))
                    sample = self.config.sample
                else:
                    sample = self.config.sampleDict[label]

            self.log.info('Creating prior for label %s, with sample %f' % (label,sample))
            if numpy.any(numpy.isnan(cov)):
                self.log.warn('Covariance matrix has nan, skipping %s'%label)
                continue

            var = float(numpy.median(exposure.getMaskedImage().getVariance().getArray()))
            ngood = 0
            sigmaFlux = numpy.sqrt(cov[0,0])
            momentPrior = lsst.meas.extensions.bfd.MomentPrior(self.config.snMin*sigmaFlux,
                                                               self.config.snMax*sigmaFlux,
                                                               cov, False,
                                                               self.config.noiseFactor,
                                                               self.config.priorSigmaCutoff,
                                                               self.config.priorSigmaStep,
                                                               self.config.nSample,
                                                               self.config.priorSigmaBuffer,
                                                               self.config.selectionOnly)
            momentPrior.setVarianceLimits(minVar, maxVar)
            for i, (src) in enumerate(sourceCat):

                if self.config.sample > 0:
                    if numpy.random.uniform(0,1) > sample:
                        continue

                bfd_control = lsst.meas.extensions.bfd.BfdKMomentControl()
                bfd_control.sigma = self.config.sigma
                bfd_control.wIndex = self.config.wIndex
                bfd_control.maxCentroid = self.config.maxXY
                bfd_control.ignorePsf = False
                bfd_control.shift = True
                bfd_control.reCentroid = self.config.reCentroid
                bfd_control.reCentroidPsf = self.config.reCentroidPsf

                try:
                    if i%1000 == 0:
                        print i
                    if src.get('bfd.flags'):
                        continue

                    pos = src.get('bfd.center')
                    priorGalaxy = lsst.meas.extensions.bfd.PriorGalaxy(bfd_control)
                    passed = priorGalaxy.addImage(src, exposure, pos, var)

                    if passed is False:
                        continue

                    weight = 1.
                    flip = True
                    momentPrior.addPriorGalaxy(priorGalaxy, self.config.maxXY, weight,
                                               flip, src.getId())
                    ngood += 1
                except MemoryError:
                    raise  # always let MemoryError propagate up, as continuing just causes more problems
                except Exception as err:
                    self.log.warn("Error measuring source %s : %s"
                                  % (i, err))

            selectionPqr = momentPrior.selectionProbability(cov)
            deselect = selectionPqr.copy()
            deselect[0] = 1 - selectionPqr[0]
            for i in range(1,6):
                deselect[i] *= -1.
            self.log.info('Used %d galaxies in prior' % ngood)
            catalog = momentPrior.getCatalog()
            metadata = lsst.daf.base.PropertyList()
            metadata.set('cov',numpy.array(cov.flatten(),dtype=float))
            metadata.set('selectionPqr',selectionPqr.astype(numpy.float))
            metadata.set('deselectPqr',deselect.astype(numpy.float))
            metadata.set('fluxMin',self.config.snMin*sigmaFlux)
            metadata.set('fluxMax',self.config.snMax*sigmaFlux)
            metadata.set('varMin',minVar)
            metadata.set('varMax',maxVar)
            metadata.set('noiseFactor',self.config.noiseFactor)
            metadata.set('priorSigmaCutoff',self.config.priorSigmaCutoff)
            metadata.set('priorSigmaStep',self.config.priorSigmaStep)
            metadata.set('priorSigmaBuffer',self.config.priorSigmaBuffer)
            metadata.set('nsample',self.config.nSample)
            metadata.set('selectionOnly',self.config.selectionOnly)
            metadata.set('invariantCovariance',self.config.invariantCovariance)
            metadata.set('maxXY', self.config.maxXY)
            metadata.set('sigma', self.config.sigma)
            metadata.set('wIndex', self.config.wIndex)
            metadata.set('centroidName', self.config.centroidName)
            metadata.set('covFile', self.config.covFile)
            metadata.set('totalWeight', momentPrior.getTotalWeight())

            catalog.getTable().setMetadata(metadata)

            self.log.info('Created %d templates' % len(catalog))
            dataRef.dataId['label'] = label

            if self.config.selectionOnly:
                dataRef.dataId['label']+='_selection'

            if self.config.noiseFactor >1:
                dataRef.dataId['label']+='_nf_%0.2f'%self.config.noiseFactor
            dataRef.put(catalog, "prior")


    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="deep_image", level="image",
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
