#
# LSST Data Management System
# Copyright 2013-2014 LSST Corporation.
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

"""
An LSST Mapper class for managing the GREAT3 simulation directory
structure.
"""

import os
import lsst.afw.geom
import lsst.daf.persistence

class DatasetDefinition(object):
    """Base class of a hierarchy used by Great3Mapper to defined different kinds of types of objects
    to persist.
    """

    __slots__ = "template", "keys", "ranges", "python", "cpp", "storage"

    def __init__(self, template, keys, ranges=None, python=None, cpp='ignored', storage=None):
        if ranges is None: ranges = {}
        self.template = template
        self.keys = keys
        self.ranges = ranges
        self.python = python
        self.cpp = cpp
        self.storage = storage

    def makeButlerLocation(self, path, dataId, suffix=""):
        """Method called by SimpleMapping to implement a map_ method."""
        return lsst.daf.persistence.ButlerLocation(self.python, self.cpp, self.storage, [path], dataId)

    def addClosures(self, cls, name, suffix=""):
        def map(mapper, dataId, write=False):
            path = os.path.join(mapper.root, self.template.format(**dataId))
            if not write:
                newPath = mapper._parentSearch(path)
                if newPath is not None:
                    path = newPath
            return self.makeButlerLocation(path, dataId, suffix="")
        setattr(cls, "map_%s%s" % (name, suffix), map)
        def query(mapper, level, format, dataId):
            current = [[]]
            for k in format:
                if k in dataId:
                    for old in current:
                        old += [dataId[k]]
                else:
                    tmp = []
                    for old in current:
                        tmp.extend(old + [v] for v in xrange(*self.ranges.get(k, mapper.ranges[k])))
                    current = tmp
            return current
        setattr(cls, "query_%s%s" % (name, suffix), query)

class ImageDatasetDefinition(DatasetDefinition):
    """DatasetDefinition type for images, with support for _sub subimage extraction."""

    def __init__(self, template, keys, ranges=None, python="lsst.afw.image.ImageF", cpp=None,
                 storage="FitsStorage"):
        if cpp is None: cpp = python.split(".")[-1]
        DatasetDefinition.__init__(self, template, keys, ranges, python, cpp, storage)

    def makeButlerLocation(self, path, dataId, suffix=""):
        """Method called by SimpleMapping to implement a map_ method; overridden to support subimages."""
        if not suffix:
            loc = DatasetDefinition.makeButlerLocation(self, path, dataId, suffix="")
        elif suffix == "_sub":
            subId = dataId.copy()
            bbox = subId.pop('bbox')
            loc = DatasetDefinition.makeButlerLocation(path, subId, suffix="")
            loc.additionalData.set('llcX', bbox.getMinX())
            loc.additionalData.set('llcY', bbox.getMinY())
            loc.additionalData.set('width', bbox.getWidth())
            loc.additionalData.set('height', bbox.getHeight())
            if 'imageOrigin' in dataId:
                loc.additionalData.set('imageOrigin',
                                       dataId['imageOrigin'])
        return loc

    def addClosures(self, cls, name, suffix=""):
        DatasetDefinition.addClosures(self, cls, name, suffix="")
        DatasetDefinition.addClosures(self, cls, name, suffix="_sub")

class CatalogDatasetDefinition(DatasetDefinition):

    def __init__(self, template, keys, ranges=None, python="lsst.afw.table.BaseCatalog", cpp=None,
                 storage="FitsCatalogStorage"):
        if cpp is None: cpp = python.split(".")[-1]
        DatasetDefinition.__init__(self, template, keys, ranges, python, cpp, storage)

class Great3Mapper(lsst.daf.persistence.Mapper):

    PIXEL_SCALE = 0.2*lsst.afw.geom.arcseconds
    PSTAMP_SIZE = 48
    FIELD_SIZE = 10*lsst.afw.geom.degrees
    TILE_SIZE = 2*lsst.afw.geom.degrees

    datasets = dict(
        image = ImageDatasetDefinition(
            template="image-{subfield:03d}-{epoch:01d}.fits",
            keys={"subfield": int, "epoch": int},
            ),
        starfield_image = ImageDatasetDefinition(
            template="starfield_image-{subfield:03d}-{epoch:01d}.fits",
            keys={"subfield": int, "epoch": int},
            ),
        deep_image = ImageDatasetDefinition(
            template="deep_image-{subfield:03d}-{epoch:01d}.fits",
            keys={"subfield": int, "epoch": int},
            ranges=dict(subfield=(0,5))
            ),
        deep_starfield_image = ImageDatasetDefinition(
            template="deep_starfield_image-{subfield:03d}-{epoch:01d}.fits",
            keys={"subfield": int, "epoch": int},
            ranges=dict(subfield=(0,5))
            ),
        galaxy_catalog = CatalogDatasetDefinition(
            template="galaxy_catalog-{subfield:03d}.fits",
            keys={"subfield": int}
            ),
        deep_galaxy_catalog = CatalogDatasetDefinition(
            template="deep_galaxy_catalog-{subfield:03d}.fits",
            keys={"subfield": int},
            ranges=dict(subfield=(0,5))
            ),
        star_catalog = CatalogDatasetDefinition(
            template="star_catalog-{subfield:03d}.fits",
            keys={"subfield": int}
            ),
        deep_star_catalog = CatalogDatasetDefinition(
            template="deep_star_catalog-{subfield:03d}.fits",
            keys={"subfield": int},
            ranges=dict(subfield=(0,5))
            ),
        src = CatalogDatasetDefinition(
            template="src-{subfield:03d}.fits",
            python="lsst.afw.table.SourceCatalog",
            keys={"subfield": int}
            ),
        deep_src = CatalogDatasetDefinition(
            template="deep_src-{subfield:03d}.fits",
            python="lsst.afw.table.SourceCatalog",
            keys={"subfield": int},
            ranges=dict(subfield=(0,5))
            ),
        shear = CatalogDatasetDefinition(
            template="shear.fits",
            python="lsst.afw.table.BaseCatalog",
            keys={}
            ),
        deep_shear = CatalogDatasetDefinition(
            template="deep_shear.fits",
            python="lsst.afw.table.BaseCatalog",
            keys={},
            ),
        star_index = CatalogDatasetDefinition(
            template="star_index-{field:03d}.fits",
            python="lsst.afw.table.BaseCatalog",
            keys={"field": int},
            )
        )

    levels = dict(
        image = [],
        subfield = ['epoch'],
        epoch = ['subfield'],
    )

    defaultLevel = "image"
    defaultSubLevels = dict(
        image=None,
        subfield="epoch"
    )

    ranges = dict(
        field = (0, 200, 20),
        subfield = (0, 200),
        epoch = (0, 1),
        )

    deep_ranges = dict(
        subfield = (0, 5),
        epoch = (0, 1),
        )

    def __init__(self, root, outputRoot=None, calibRoot=None):  # calibRoot ignored
        # this block copied directly from CameraMapper.__init__ in daf_butlerUtils
        if outputRoot is not None:
            # Path manipulations are subject to race condition
            if not os.path.exists(outputRoot):
                try:
                    os.makedirs(outputRoot)
                except OSError, e:
                    if not e.errno == errno.EEXIST:
                        raise
                if not os.path.exists(outputRoot):
                    raise RuntimeError, "Unable to create output " \
                            "repository '%s'" % (outputRoot,)
            if os.path.exists(root):
                # Symlink existing input root to "_parent" in outputRoot.
                src = os.path.abspath(root)
                dst = os.path.join(outputRoot, "_parent")
                if not os.path.exists(dst):
                    try:
                        os.symlink(src, dst)
                    except OSError:
                        pass
                if os.path.exists(dst):
                    if os.path.realpath(dst) != os.path.realpath(src):
                        raise RuntimeError, "Output repository path " \
                                "'%s' already exists and differs from " \
                                "input repository path '%s'" % (dst, src)
                else:
                    raise RuntimeError, "Unable to symlink from input " \
                            "repository path '%s' to output repository " \
                            "path '%s'" % (src, dst)
            # We now use the outputRoot as the main root with access to the
            # input via "_parent".
            root = outputRoot

        self.root = root
        self.index = {}

    def getDefaultLevel(self): return self.defaultLevel

    def getDefaultSubLevel(self, level):
        return self.defaultSubLevels.get(level, None)

    def getKeys(self, datasetType, level):
        if datasetType is None:
            keyDict = {}
        else:
            keyDict = self.datasets[datasetType].keys
        if level is not None and level in self.levels:
            keyDict  = dict(keyDict)
            for l in self.levels[level]:
                if l in keyDict:
                    del keyDict[l]
        return keyDict

    # More stuff copied directly from CameraMapper
    def _parentSearch(self, path):

        # Separate path into a root-equivalent prefix (in dir) and the rest
        # (left in path)
        rootDir = self.root

        # First remove trailing slashes (#2527)
        while len(rootDir) > 1 and rootDir[-1] == '/':
            rootDir = rootDir[:-1]

        if path.startswith(rootDir + "/"):
            # Common case; we have the same root prefix string
            path = path[len(rootDir)+1:]
            dir = rootDir
        elif rootDir == "/" and path.startswith("/"):
            path = path[1:]
            dir = rootDir
        else:
            # Search for prefix that is the same as root
            pathPrefix = os.path.dirname(path)
            while pathPrefix != "" and pathPrefix != "/":
                if os.path.realpath(pathPrefix) == os.path.realpath(self.root):
                    break
                pathPrefix = os.path.dirname(pathPrefix)
            if os.path.realpath(pathPrefix) != os.path.realpath(self.root):
                # No prefix matching root, don't search for parents
                if os.path.exists(path):
                    return path
                return None
            if pathPrefix == "/":
                path = path[1:]
            elif pathPrefix != "":
                path = path[len(pathPrefix)+1:]
            # If pathPrefix == "", then the current directory is the root
            dir = pathPrefix

        # Now search for the path in the root or its parents
        while not os.path.exists(os.path.join(dir, path)):
            dir = os.path.join(dir, "_parent")
            if not os.path.exists(dir):
                return None
        return os.path.join(dir, path)

    # More stuff copied directly from CameraMapper
    def backup(self, datasetType, dataId):
        n = 0
        suffix = ""
        newLocation = self.map(datasetType, dataId, write=True)
        newPath = newLocation.getLocations()[0]
        path = self._parentSearch(newPath)
        oldPaths = []
        while path is not None:
            n += 1
            oldPaths.append((n, path))
            path = self._parentSearch("%s~%d" % (newPath, n))
        for n, oldPath in reversed(oldPaths):
            newDir, newFile = os.path.split(newPath)
            if not os.path.exists(newDir):
                os.makedirs(newDir)
            shutil.copy(oldPath, "%s~%d" % (newPath, n))

    @staticmethod
    def getCameraName():
        return "great3"

    @staticmethod
    def getEupsProductName():
        return "obs_great3"

def addClosures(cls):
    for name, dataset in cls.datasets.iteritems():
        dataset.addClosures(cls, name)

addClosures(Great3Mapper)
