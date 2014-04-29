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

class FieldTaskRunner(lsst.pipe.base.TaskRunner):
    """A TaskRunner that passes DataRefs to the run() method by tile.
    """

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        byField = {}
        for dataId in parsedCmd.id.idList:
            dataRef = parsedCmd.butler.dataRef("star_catalog", **dataId)
            field = (dataId['subfield'] // 20) * 20
            byField.setdefault(field , []).append(dataRef)
        return [(parsedCmd.butler.dataRef("star_index", field=field), dataRefList, kwargs)
                for field, dataRefList in byField.iteritems()]

    def __call__(self, args):
        fieldDataRef, dataRefList, kwargs = args
        task = self.makeTask(args=args)
        task.run(fieldDataRef, dataRefList, **kwargs)


class MakeStarIndexConfig(lsst.pex.config.Config):
    pass

class MakeStarIndexTask(lsst.pipe.base.CmdLineTask):

    RunnerClass = FieldTaskRunner
    ConfigClass = MakeStarIndexConfig

    _DefaultName = "makeStarIndex"

    def __init__(self, **kwds):
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        self.schema = lsst.afw.table.Schema()
        self.subfieldKey = self.schema.addField("subfield.index", type=long, doc="subfield index")
        self.tileIndexKeyX = self.schema.addField("tile.index.x", type=long, doc="tile index x")
        self.tileIndexKeyY = self.schema.addField("tile.index.y", type=long, doc="tile index y")
        self.subfieldPosKeyX = self.schema.addField("subfield.position.x", type=long,
                                                    doc="x position on original subfield grid image")
        self.subfieldPosKeyY = self.schema.addField("subfield.position.y", type=long,
                                                    doc="y position on original subfield grid image")
        self.tilePosKeyX = self.schema.addField("tile.position.x", type=long,
                                                doc="x position on final tile image")
        self.tilePosKeyY = self.schema.addField("tile.position.y", type=long,
                                                doc="y position on final tile image")

    def run(self, fieldDataRef, dataRefList):
        PIXEL_SCALE = fieldDataRef.butlerSubset.butler.mapper.PIXEL_SCALE
        index = lsst.afw.table.BaseCatalog(self.schema)
        mapper = None
        for dataRef in dataRefList:
            self.log.info("Processing stars from subfield=%s" % dataRef.dataId['subfield'])
            cat = dataRef.get("star_catalog", immediate=True)
            if mapper is None:
                mapper = lsst.afw.table.SchemaMapper(cat.schema, self.schema)
                mapper.addMapping(cat.schema.find("x.tile.index").key, "tile.index.x", True)
                mapper.addMapping(cat.schema.find("y.tile.index").key, "tile.index.y", True)
                mapper.addMapping(cat.schema.find("x").key, "subfield.position.x", True)
                mapper.addMapping(cat.schema.find("y").key, "subfield.position.y", True)
                xTilePosKey = cat.schema.find("tile.x.pos.deg").key
                yTilePosKey = cat.schema.find("tile.y.pos.deg").key
            tmp = lsst.afw.table.BaseCatalog(index.table)
            tmp.extend(cat, mapper=mapper)
            tmp[self.subfieldKey][:] = dataRef.dataId['subfield']
            # The true degree positions of the stars from GREAT3 don't put things on the
            # image grid for the tile, so we're forced to round off the positions of the
            # stars.  Hopefully, this won't matter, because the PSF shouldn't be varying
            # on subpixel scales.  Alternatively, we could resample, but that's slower
            # and perhaps more dangerous in terms of introducing systematics.
            numpy.rint(cat[xTilePosKey] / PIXEL_SCALE.asDegrees(), out=tmp[self.tilePosKeyX])
            numpy.rint(cat[yTilePosKey] / PIXEL_SCALE.asDegrees(), out=tmp[self.tilePosKeyY])
            index.extend(tmp, deep=False, mapper=None)
        index.sort(self.tileIndexKeyX)
        index.sort(self.tileIndexKeyY)
        fieldDataRef.put(index, "star_index")

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument(name="--id", datasetType="star_catalog", level="image",
                               help="data ID, e.g. --id subfield=0..20")
        return parser

    def writeConfig(self, butler, clobber=False):
        pass

    def writeSchemas(self, butler, clobber=False):
        pass

    def writeMetadata(self, dataRef):
        pass

    def writeEupsVersions(self, butler, clobber=False):
        pass
