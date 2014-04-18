try:
    import lsst.meas.extensions.shapeHSM
    root.measurement.algorithms.names |= ["shape.hsm.regauss"]
except ImportError:
    print "shapeHSM not found; running only basic algorithms"
