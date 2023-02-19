import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import plotNcFunc as pf

import multiprocessing
from multiprocessing import Pool

leap = 400
to = 100000
if __name__ == '__main__':  
    try:
        nProc = int(multiprocessing.cpu_count() / 2)
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotOnCubeWindMul, (t, )) for t in range(0, to, leap)]
            final = [result.get() for result in results]
    except:
        print("finish1")

    try:
        nProc = int(multiprocessing.cpu_count() / 2)
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotOnSphereWindMul, (t, )) for t in range(0, to, leap)]
            final = [result.get() for result in results]
    except:
        print("finish2")

    try:
        nProc = int(multiprocessing.cpu_count() / 2)
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotSphereWindCartopy, (t, )) for t in range(0, to, leap)]
            final = [result.get() for result in results]
    except:
        print("finish3")
