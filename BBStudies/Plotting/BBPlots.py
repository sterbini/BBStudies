

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import pandas as pd
import scipy.special as sciSpec
import numpy as np
import itertools
from fractions import Fraction

import BBStudies.Physics.Resonances as resn

#############################################################
# WORKING DIAGRAM

def _isPointInside(x,y,x_range,y_range,tol = 0):
    x1,x2 = x_range
    y1,y2 = y_range
    tol   = np.abs(tol)
    return (x1-tol <= x <= x2+tol) and (y1-tol <= y <= y2+tol)

def _isLineInside(a,b,x_range,y_range,tol = 0):
    '''Assumes y = ax+b and finds if line goes inside the ROI'''
    x1,x2 = np.sort(x_range)
    y1,y2 = np.sort(y_range)
    tol   = np.abs(tol)

    xi1   = (y1-b)/a if a!=0 else 0
    xi2   = (y2-b)/a if a!=0 else 0

    yi1   = a*x1 + b
    yi2   = a*x2 + b

    return ((x1-tol <= xi1 <= x2+tol) or 
            (x1-tol <= xi2 <= x2+tol) or 
            (y1-tol <= yi1 <= y2+tol) or 
            (y1-tol <= yi2 <= y2+tol))

def _plot_resonance_lines(df,Qx_range,Qy_range,ROI_tol = 1e-3,**kwargs):
    ''' Plots all resonance lines contained in a dataframe, provided by BBStudies.Physics.Resonances.resonance_df
        -> Filters out the lines that do not enter the ROI defined by Qx_range,Qy_range'''

    options = {'color':'k','alpha':0.15}
    options.update(kwargs)

    if 'label' in options.keys():
        _label = options.pop('label')
        plt.plot([np.nan],[np.nan],label = _label,**options)

    for _,line in df.iterrows():
        
        # Non-vertical lines
        if line['slope'] != np.inf:
            # Skip if line not in ROI
            if not _isLineInside(line['slope'],line['y0'],Qx_range,Qy_range,tol = ROI_tol):
                continue

            xVec = np.array(Qx_range)
            yVec = line['slope']*xVec + line['y0']

        # Vertical line
        else:
            # Skip if line not in ROI
            if not (Qx_range[0] <= line['x0'] <=Qx_range[1]):
                continue

            xVec = line['x0']*np.ones(2)
            yVec = np.array(Qy_range)

        # Plotting if in ROI
        plt.plot(xVec,yVec,**options)





def workingDiagram(Qx_range = [0,1],Qy_range = [0,1],order=6,**kwargs):
    
    if not isinstance(order, (list, type(np.array([])))):
        # Regular, full working diagram
        resonances = resn.resonance_df(order)
        _plot_resonance_lines(resonances,Qx_range,Qy_range,ROI_tol = 1e-3,**kwargs)
    else:
        # Selected resonances
        resonances = resn.resonance_df(np.max(order))
        for _ord in order:
            _plot_resonance_lines(resonances[resonances['Order']==_ord],Qx_range,Qy_range,ROI_tol = 1e-3,**kwargs)

        


   


    
   


def xticks_from_farey(order,*args,**kwargs):
    plt.xticks([m/k for m,k in resn.Farey(order)],[f'{m:d}/{k:d}' for m,k in resn.Farey(order)],*args,**kwargs)


def yticks_from_farey(order,*args,**kwargs):
    plt.yticks([m/k for m,k in resn.Farey(order)],[f'{m:d}/{k:d}' for m,k in resn.Farey(order)],*args,**kwargs)
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


#############################################################
def FFT(x,*args,flipped=False,unpack=False,**kwargs):
    x     = np.array(x)
    turns = np.arange(1,len(x)+1)
    freq  = np.fft.fftfreq(turns.shape[-1])

    spectrum = np.fft.fft(x-np.mean(x))
    #idx      = np.argmax(np.abs(spectrum))
    #Qx       = freq[idx]
    if flipped:
        plt.plot(freq[freq>0],-np.abs(spectrum)[freq>0],*args,**kwargs)
    else:
        plt.plot(freq[freq>0],np.abs(spectrum)[freq>0],*args,**kwargs)

    if unpack:
        return freq[freq>0],np.abs(spectrum)[freq>0]
#############################################################