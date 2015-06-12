try:
    import lsst.meas.extensions.shapeHSM
    root.measurement.algorithms.names.add("ext_shapeHSM_HsmShapeRegauss")
except ImportError:
    print "shapeHSM not found; running only basic algorithms"
