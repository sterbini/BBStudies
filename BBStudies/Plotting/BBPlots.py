

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import pandas as pd
import scipy.special as sciSpec
import numpy as np
import itertools


#############################################################
def plotWorkingDiagram(order = 12,QxRange=np.array([0,1]),QyRange=np.array([0,1]),**kwargs):

    # Initialization
    options = {'color':'k','alpha':0.5}
    options.update(kwargs)
    QxRange,QyRange = np.array(QxRange),np.array(QyRange)
    def intList(n): return np.arange(-n,n+1)
    plt.axis('square')
    plt.xlim(QxRange)
    plt.ylim(QyRange)


    # Creating all combinations except vertical lines
    popt = []
    for m1,m2,n in itertools.product(intList(order),intList(order)[intList(order)!=0],intList(200)):
        if np.abs(m1)+np.abs(m2) <= order:
            popt.append((-m1/m2,n/m2))


    # Removing Duplicates
    # TODO: change this line in order to keep track of the order of the resonance
    popt = list(set(popt))

    # Keeping only lines in ROI
    ROI_popt = []
    for slope,y0 in popt:
        line = slope*QxRange + y0

        if np.any(np.logical_and(line>=np.min(QyRange),line<=np.max(QyRange))):
            ROI_popt.append((slope,y0))

    # Plotting
    regularSlopes = np.array(ROI_popt)[:,0]
    for slope,y0 in ROI_popt:
        plt.plot(QxRange,slope*QxRange + y0,**options)

        # Reflection around y=x to take care of the cases where m2=0
        with np.errstate(divide='ignore'):
            invertedSlope = (np.diff(QyRange)/np.diff(slope*QyRange + y0))[0]
        if not np.round(invertedSlope,5) in list(np.round(regularSlopes,5)):
            plt.plot(slope*QyRange + y0,QyRange,**options)

#############################################################


#############################################################
def boundedScatter(x,y,c,boundaries,cmap='viridis',zorder=2, **kwargs):
    norm = BoundaryNorm(boundaries= boundaries, ncolors=int(0.9*256))
    sc   = plt.scatter(x,y,c=c,norm=norm,zorder=2,**kwargs)

    return sc
#############################################################



#############################################################
def polarmesh(x,y,r,theta,*args,**kwargs):
    _df = pd.DataFrame({'x':np.array(x),'y':np.array(y),'r':np.array(r),'theta':np.array(theta)})

    options = {'color':'darkslateblue','alpha':0.3}
    options.update(kwargs)

    for sortKey in ['r','theta']:
        for name, group in _df.groupby(sortKey):
            plt.plot(group['x'],group['y'],*args,**options)

#############################################################