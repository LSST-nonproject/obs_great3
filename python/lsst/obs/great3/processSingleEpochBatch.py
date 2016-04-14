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
import lsst.pipe.base
from hsc.pipe.base.pool import abortOnError, Pool
from hsc.pipe.base.parallel import BatchPoolTask
import hsc.pipe.base.butler as hscButler

from .processSingleEpoch import *

class Runner(lsst.pipe.base.TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """Get bare butler into Task"""
        if len(parsedCmd.id.refList) > 0 and len(parsedCmd.deep_id.refList) > 0:
            raise RuntimeError("You can only specify id or deep_id, not both")
        if len(parsedCmd.id.refList) > 0:
               refList = parsedCmd.id.refList
        else:
               refList = parsedCmd.deep_id.refList
        return [(refList, kwargs),]

class ProcessSingleEpochBatchConfig(lsst.pex.config.Config):
    singleEpoch = lsst.pex.config.ConfigurableField(
        target = ProcessSingleEpochTask,
        doc = "Single epoch processing task to run on multiple cores"
    )

class ProcessSingleEpochBatchTask(BatchPoolTask):

    ConfigClass = ProcessSingleEpochBatchConfig
    RunnerClass = Runner
    _DefaultName = "processSingleEpochBatch"

    def __init__(self, *args, **kwargs):
        super(ProcessSingleEpochBatchTask, self).__init__(*args, **kwargs)
        self.makeSubtask("singleEpoch")

    @abortOnError
    def run(self, dataRefList):
        """Process a single exposure, with scatter-gather-scatter using MPI.

        All nodes execute this method, though the master and slaves have different
        routes through it.  The expRef is only a DummyDataRef on the slaves.
        """
        pool = Pool("processSingleEpochBatch")
        pool.cacheClear()
        resultList = pool.map(self.process, dataRefList)

    def process(self, cache, dataRef):
        with self.logOperation("processing %s" % (dataRef.dataId)):
            try:
                result = self.singleEpoch.run(dataRef)
            except Exception, e:
                self.log.warn("Failed to process %s: %s\n" % (dataRef.dataId, e))
                import traceback
                traceback.print_exc()
                return None
        return result

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numNodes, numProcs):
        return time

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        doBatch = kwargs.pop("doBatch", False)
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="image", level="image",
                               help="data ID, e.g. --id subfield=0")
        parser.add_id_argument(name="--deep_id", datasetType="deep_image", level="image",
                               help="deep_data ID, e.g. --deep_id subfield=0")

        return parser

    def writeConfig(self, butler, clobber=False, doBackup=False):
        pass

    def writeSchemas(self, butler, clobber=False, doBackup=False):
        pass

    def writeMetadata(self, dataRef):
        pass

    def writeEupsVersions(self, butler, clobber=False, doBackup=False):
        pass

    def _getConfigName(self):
        return None

    def _getEupsVersionsName(self):
        return None

    def _getMetadataName(self):
        return None
