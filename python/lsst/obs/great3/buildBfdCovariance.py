#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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



class Runner(lsst.pipe.base.TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """Get bare butler into Task"""
        return [(parsedCmd.id.refList, kwargs),]

class BuildBfdCovarianceConfig(lsst.pex.config.Config):
     bins = lsst.pex.config.Field(
        dtype=int,
        default=5,
        optional=True,
        doc="Split data set into this many bins and compute a covariance matrix for each one"
     )
     label = lsst.pex.config.Field(
        dtype=str,
        default='b',
         optional=True,
         doc="This will be the labe for each bin"
     )
     minPct = lsst.pex.config.Field(
         dtype=float,
         default=0.1,
         optional=True,
         doc="The lower bound of noise variance to consider"
     )
     maxPct = lsst.pex.config.Field(
         dtype=float,
         default=99.9,
         optional=True,
         doc="The upper bound of noise variance to consider"
     )
     outputFile = lsst.pex.config.Field(
         dtype=str,
         default='./covariance.fits',
         optional=True,
         doc="Where to store the output"
     )

class BuildBfdCovarianceTask(lsst.pipe.base.CmdLineTask):

    ConfigClass = BuildBfdCovarianceConfig
    _DefaultName = "buildBfdCovariance"
    RunnerClass = Runner


    def __init__(self,  **kwargs):
        """Initialize
        """
        lsst.pipe.base.CmdLineTask.__init__(self, **kwargs)

        self.schema = lsst.afw.table.Schema()
        self.labelKey = self.schema.addField("label",type=str, doc="name of bin",size=10)
        self.minKey = self.schema.addField("min",type=float, doc="minimum value of the variance")
        self.maxKey = self.schema.addField("max",type=float, doc="maximum value of the variance")
        self.isoCovKey = self.schema.addField("isoCov", doc="isotropized moment covariance matrix",
                                           type="ArrayF", size=36)
        self.covKey = self.schema.addField("cov", doc="moment covariance matrix", type="ArrayF", size=36)
        self.catalog =  lsst.afw.table.BaseCatalog(self.schema)


    def run(self, dataRefList):

        noiseVariance=[]
        cov =[]
        for dataRef in dataRefList:
            self.log.info('Adding %s'%str(dataRef.dataId))
            try:
                cat = dataRef.get('src', flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)
                maskGood = cat.get('bfd.flags') == False
                noiseVariance.append(cat.get('bfd.momentsCov')[:,0][maskGood][0])
                cov.append(cat.get('bfd.momentsCov')[maskGood][0])
            except Exception as e:
                self.log.info('Failed for %s: %s'% (str(dataRef.dataId),e))

        cov=numpy.array(cov)
        noiseVariance=numpy.array(noiseVariance)
        minVariance = numpy.percentile(noiseVariance, self.config.minPct)
        maxVariance = numpy.percentile(noiseVariance, self.config.maxPct)
        print noiseVariance
        self.log.info('Min/max variances %0.4f/%0.4f'% (minVariance,maxVariance))

        binWidth = (maxVariance - minVariance) / self.config.bins
        bins = numpy.arange(minVariance, maxVariance + binWidth, binWidth)

        for i in range(len(bins)-1):
            rec = self.catalog.addNew()
            mask = numpy.logical_and(noiseVariance >= bins[i], noiseVariance <= bins[i+1])

            tot = numpy.zeros(36)
            for val in cov[mask]:
                tot += val
            tot /= numpy.sum(mask)
            rec.set('cov',numpy.array(tot.flatten(),dtype=numpy.float32))
            #
            # isotropize matrix
            sum = tot.reshape(6,6)

            varE = 0.5*(sum[4,4] + sum[5,5])
            varC = 0.5*(sum[1,1] + sum[2,2])
            sum[1:3,:]=0
            sum[:,1:3]=0
            sum[4:6,:]=0
            sum[:,4:6]=0
            sum[1,1]=varC
            sum[2,2]=varC
            sum[4,4]=varE
            sum[5,5]=varE
            print bins[i],bins[i+1],numpy.sum(mask)
            rec.set('max',bins[i+1])
            rec.set('min',bins[i])
            rec.set('isoCov',numpy.array(sum.flatten(),dtype=numpy.float32))

            rec.set('label', '%s%d' % (self.config.label, i))

        self.catalog.writeFits(self.config.outputFile)



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

