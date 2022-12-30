import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import plotFunc as pf

import multiprocessing
from multiprocessing import Pool

DX = DY = 2
NX = NY = int(90 / DX + 2- 2)
DT = 360

if __name__ == '__main__':
    try:
        nProc = multiprocessing.cpu_count() - 4
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotOnCubeWindMul, (t, )) for t in range(0, 1000, 10)]
            final = [result.get() for result in results]
    except:
        print("finish1")

    try:
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotOnSphereWindMul, (t, )) for t in range(0, 1000, 10)]
            final = [result.get() for result in results]
    except:
        print("finish2")

    try:
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotSphereWindCartopy, (t, )) for t in range(0, 1000, 10)]
            final = [result.get() for result in results]
    except:
        print("finish3")