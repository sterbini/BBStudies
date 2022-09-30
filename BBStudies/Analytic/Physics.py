import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools


import BBStudies.Base.Constants as cst




##############################################################
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)
##############################################################


##############################################################
def vecNorm(vec):
    return np.sqrt(vec[0]**2 + vec[1]**2)
##############################################################


##############################################################
def normCoordinates(x,px,alpha=None,beta=None,SVD=False):
    
    if SVD:
        alpha,beta = SVD_AlphaBeta(x,px)
        
    #N0 = [[1/np.sqrt(beta),0],[alpha/np.sqrt(beta), np.sqrt(beta)]]
    x_n  = x/np.sqrt(beta)
    px_n = alpha*x/np.sqrt(beta) + px*np.sqrt(beta)
    
    return x_n,px_n
##############################################################


##############################################################
def SVD_AlphaBeta(x,px):
    '''Taken from https://arxiv.org/pdf/2006.10661.pdf '''
    
    U,s,V= np.linalg.svd([x,px])         #SVD
    
    N = np.dot(U,np.diag(s))
    theta = np.arctan(-N[0,1]/N[0,0])    #AngleofR(theta)
    co=np.cos(theta) ; si=np.sin(theta)
    
    R = [[co,si],[-si,co]]   
    X = np.dot(N,R)                      #Floquetupto1/det(USR)
    
    beta = np.abs(X[0,0]/X[1,1])
    alpha = X[1,0]/X[1,1]
    
    # dropped
    ex =s[0]*s[1]/(len(x)/2.)            #emit=det(S)/(n/2)
    
    return alpha,beta

##############################################################


##############################################################
def getAction(x,px,alpha=None,beta=None,SVD=False):
    
    if SVD:
        alpha,beta = computeAlphaBeta(x,px)
    gamma = (1+alpha**2)/beta
    
    J = (gamma*x**2  + 2*alpha*x*px + beta*px**2)/2
    
    return J
##############################################################


##############################################################
def generate_coordGrid(xRange,yRange,labels = ['x','y'],nPoints=100):
    '''
    Distribute points uniformly on a 2D grid.
    -----------------------------------------
    Input:
        xRange : range of first coordinate
        yRange : range of second coordinate
        labels : labels to be used in the resulting dataframe
        nPoint : total number of points to generate (sqrt(nPoints) for each coordinate)
    Returns:
        coordinates: dataframe containing the distributed points
    '''

    if type(xRange) is list and type(yRange) is list:
        xVec = np.linspace(xRange[0],xRange[1],int(np.sqrt(nPoints)))
        yVec = np.linspace(yRange[0],yRange[1],int(np.sqrt(nPoints)))
    else:
        xVec = xRange
        yVec = yRange
        
    xx,yy = np.meshgrid(xVec,yVec)
    xx,yy = xx.flatten(),yy.flatten()

    return pd.DataFrame(dict(zip(labels,[xx,yy])))
##############################################################


##############################################################
def polar_grid(r_n,theta_n,emitt = None):
    _coord  = generate_coordGrid(r_n,theta_n,labels = ['r_n','theta_n'])

    _coord.insert(2,'x_n',_coord['r_n']*np.cos(_coord['theta_n']))
    _coord.insert(3,'y_n',_coord['r_n']*np.sin(_coord['theta_n']))

    if emitt is not None:
        if isinstance(emitt, (int, float)):
            emitt = [emitt,emitt]
        _coord.insert(4,'J_x',(_coord['x_n']**2)*emitt[0]/2)
        _coord.insert(5,'J_y',(_coord['y_n']**2)*emitt[1]/2)

    return _coord
##############################################################