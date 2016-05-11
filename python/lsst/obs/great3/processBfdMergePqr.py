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

import lsst.pex.config
import lsst.pipe.base
import lsst.afw.image
import lsst.afw.geom
import lsst.afw.table
import lsst.meas.extensions.bfd
import numpy
from lsst.daf.persistence import Butler
import lsst.meas.extensions.bfd as bfd
import glob
import os
from lsst.pipe.base import Struct, CmdLineTask, ArgumentParser, TaskRunner, TaskError

from .processSingleEpoch import *

class Runner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        refList = {} # Will index this as refList[subfield][label] = ref
        for ref in parsedCmd.id.refList:
            subfield = ref.dataId["subfield"]
            label = ref.dataId["label"]

            if not subfield in refList:
                refList[subfield] = {}
            if label in refList[subfield]:
                raise RuntimeError("Multiple versions of %s" % (ref.dataId,))
            if not label in refList[subfield]:
                refList[subfield][label] = {}


            refList[subfield][label] = ref
        return [(t.values(), kwargs) for t in refList.itervalues()]

class ProcessBfdMergePqrConfig(ProcessSingleEpochConfig):
    selectFile = lsst.pex.config.Field(dtype=str, default="/tigress/rea3/hsc_sim/priorSelection.fits", doc="Name of coadd")
    snBins = lsst.pex.config.Field(dtype=int, default=2, doc="Number of S/N bins")

class ProcessBfdMergePqrTask(ProcessSingleEpochTask):

    ConfigClass = ProcessBfdMergePqrConfig
    _DefaultName = "processBfdMergePqr"
    RunnerClass = Runner


    def __init__(self,  **kwargs):
        """Initialize the measurement task, including the modelfits catalog schema,
        the model, prior, and calib objects, and the fitter subtask.
        """

        ProcessSingleEpochTask.__init__(self, **kwargs)
        selectFile = lsst.afw.table.BaseCatalog.readFits(self.config.selectFile)
        self.select = {}
        self.deselect = {}
        for rec in selectFile:
            self.select[rec.get('label')] = rec.get('selection')
            self.deselect[rec.get('label')] = rec.get('deselection')

    def run(self, dataRefList):
        """Main driver
        """
        catalogs = dict(self.readCatalog(dataRef) for dataRef in dataRefList)
        mergedCatalog = self.mergeCatalogs(catalogs)
        dataRefList[0].put(mergedCatalog, 'merge_pqr')

    def readCatalog(self, dataRef):
        """Read input catalog

        We read the input dataset provided by the 'inputDataset'
        class variable.
        """
        labelName = dataRef.dataId["label"]
        catalog = dataRef.get("pqr", immediate=True)
        self.log.info("Read %d pqr sources for label %s: %s" % (len(catalog), labelName, dataRef.dataId))
        return labelName, catalog

    def mergeCatalogs(self, catalogs):
        """Merge multiple catalogs

        catalogs: dict mapping label name to source catalog

        Returns: merged catalog
        """
        cat = catalogs[catalogs.keys()[0]]
        schema = lsst.afw.table.Schema()

        mapper = lsst.afw.table.SchemaMapper(cat.schema)
        mapper.addMinimalSchema(lsst.afw.table.SourceTable.makeMinimalSchema(),True)
        for item in cat.schema:
            mapper.addMapping(item.key, item.field.getName())

        outSchema = mapper.getOutputSchema()
        labelKey = outSchema.addField('prior_label',type=str,doc='label of prior',size=15)
        nsKey = {}
        nsPqrKey = {}
        if self.config.snBins > 1:
            for bin in range(self.config.snBins):
                nsPqrKey[bin] = outSchema.addField('bfd.ns.pqr.%d'%bin,type='ArrayF',
                                                   doc='non selection term for S/N bin %d'%bin, size=6)
                nsKey[bin] = outSchema.addField('bfd.ns.flux.%d'%bin,type='Flag',
                                                doc='non selection S/N flag bin %d'%bin)

        catalog = lsst.afw.table.SourceCatalog(outSchema)

        flagKey = cat.schema.find('bfd.flags').key
        fluxKey = cat.schema.find('bfd.ns.flux').key
        pqrKey = cat.schema.find('bfd.pqr').key
        labelNames = catalogs.keys()

        self.log.info('combining %s'% str(labelNames))

        for i,records in enumerate(zip(*catalogs.values())):

            new_record = catalog.addNew()

            # If we only have one bin then we can just copy the whole record
            if self.config.snBins==1:
                new_record.assign(records[0],mapper)
                continue

            new_record.set(flagKey,True)
            if self.config.snBins > 1:
                for bin in range(self.config.snBins):
                    new_record.set(nsKey[bin],False)

            goodMeasure = False
            # get list of all measurements that were within the flux limits and had no flags
            flags = numpy.array([((record.get(flagKey) == False) & (record.get(fluxKey) == False))
                                 for record in records])

            if numpy.sum(flags) > 1:
                self.log.warn("Found more than one measurement for this record, this shouldn't happen")
                continue

            index = numpy.where(flags==True)[0]
            if len(index) == 1:
                new_record.set(labelKey,labelNames[index[0]])
                new_record.assign(records[index[0]],mapper)
                new_record.set(flagKey, False)

            # Identify objects that are outside of flux limits
            ns_flags = numpy.array([ ((record.get(flagKey) == False) & (record.get(fluxKey) == True))
                                     for record in records])
            index = numpy.where(ns_flags==True)[0]

            # if only one record is outside of the flux limit assign selection Pqr
            if numpy.sum(ns_flags) == 1:
                name = labelNames[index[0]]
                if len(name)<4:
                    snBin=0
                else:
                    snBin=1
                new_record.set(nsPqrKey[snBin], records[index[0]].get(pqrKey))
                new_record.set(nsKey[snBin], True)

            # Identify objects that are not in the flux selection region of all S/N bins
            if numpy.sum(ns_flags) > 1:
                names = numpy.array(labelNames)[ns_flags]
                lnames = [len(a) for a in names]
                ns_label = numpy.array(names)[lnames == numpy.min(lnames)][0]
                new_record.set(pqrKey, self.deselect[ns_label])
                new_record.set(fluxKey, True)
                new_record.set(flagKey, False)
                new_record.set(labelKey,ns_label)

        return catalog

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="pqr", level="image",
                               help="data ID, e.g. --id subfield=0 label=b1^b2")
        return parser


