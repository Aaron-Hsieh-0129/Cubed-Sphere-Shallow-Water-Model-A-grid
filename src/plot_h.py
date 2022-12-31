import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import plotFunc as pf

import multiprocessing
from multiprocessing import Pool

leap = 20
if __name__ == '__main__':
    pf.plotWind()
    try:
        nProc = multiprocessing.cpu_count() - 4
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotOnCubeMul, (t, )) for t in range(0, 10000, leap)]
            final = [result.get() for result in results]
    except:
        print("finish1")

    try:
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotOnSphereMul, (t, )) for t in range(0, 10000, leap)]
            final = [result.get() for result in results]
    except:
        print("finish2")

    try:
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotSphereCartopy, (t, )) for t in range(0, 10000, leap)]
            final = [result.get() for result in results]
    except:
        print("finish3")
