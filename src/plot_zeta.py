import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import plotFunc as pf

import multiprocessing
from multiprocessing import Pool

leap = 20
if __name__ == '__main__':
    try:
        nProc = multiprocessing.cpu_count() - 4
        with Pool(nProc) as p:
            results = [p.apply_async(pf.plotSphereCartopyZeta, (t, )) for t in range(0, 10000, leap)]
            final = [result.get() for result in results]
    except:
        print("finish1")
