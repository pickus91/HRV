# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 19:19:25 2017

@author: Sarah Pickus

Poincare plots are an important visualization technique for quantifying the non-linear
characteristics of the RR interval time series. They are generated via plotting 
each RR interval (RR[n]) against the subsequent RR interval (RR[n+1]). There are 
a number of geometrical descriptors that can be extracted from the poincare
plot that have been shown to provide insights into short and long term HRV
trends. This package allows the user to calculate a number of these poincare
features, which are detailed below:

    elipseFittingMethod   : The ellipse fitting technique is one of the most popular
                            methods for characterizing the geometry of the poincare 
                            plot and providing insights into the level of short- and 
                            long-term variability. The short-term variability is
                            represented by SD1, which is the standard deviation of points 
                            perpendicular to the axis line-of-identy, while the 
                            long-term variability is represented by SD2, which is the 
                            tandard deviation of points along the line-of-identy.                           
                            This method is described in detail in:
                                
                            Brennan, Michael, Marimuthu Palaniswami, and Peter Kamen. 
                            "Do existing measures of Poincare plot geometry reflect 
                            nonlinear features of heart rate variability?." IEEE 
                            transactions on biomedical engineering 48.11 (2001): 
                            1342-1347.
        
    hraMethod             : C_UP and C_DOWN are derived from the standard descriptor
                            SD1^2, the variance of which corresponds to short-term HRV.
                            These more refined descriptors quantify the Poincare plot 
                            assymetry about the line-of-identity and in doing so 
                            characterize accelerations (C_DOWN) and decelerations (C_UP) 
                            of heart rate. More details on this method can be found here:
                            
                            Piskorski, J., and P. Guzik. "Geometry of the Poincaré plot 
                            of RR intervals and its asymmetry in healthy adults." 
                            Physiological measurement 28.3 (2007): 287.
        
    correlationCoefficient: The interbeat autocorrelation coefficient of the RR
                            interval time series (Pearson's correlation coefficient) 
                            has been shown to characterize poincare plot geometry,
                            as shown in:
                                
                            Otzenberger, Hélène, et al. "Dynamic heart rate variability: 
                            a tool for exploring sympathovagal balance continuously 
                            during sleep in men." American Journal of Physiology 
                            275 (1998): H946-H950.

Additionally, this package allows you to generate a poincare plot and generate histogram
distributions from mulitple projections. Statistical properties of these projections
can then be calculated as additional descriptors. 

"""

import numpy as np
from matplotlib import style
import matplotlib.pyplot as plt
style.use('ggplot')

def plotPoincare(RRints):
    """
    Input    :
    
     - RRints: [list] of RR intervals
        
    Output   :

     - Poincare plot     
    """
    ax1 = RRints[:-1]
    ax2 = RRints[1:]   
    plt.scatter(ax1, ax2, c = 'r', s = 12)
    plt.xlabel('RR_n (s)')
    plt.ylabel('RR_n+1 (s)')
    plt.show()

def eclipseFittingMethod(RRints):
    """
    Input        :
    
     - RRints    : [list] of RR intervals
        
    Output       : 
             
      - SD1, SD2 : {dict} with keys 'SD1' (numpy.float64), representing short-term 
                   variation, and 'SD2' (numpy.float64), representing long-term
                   variation.   
    """
    SDSD =  np.std(np.diff(RRints))
    SDRR = np.std(RRints)
    SD1 = (1 / np.sqrt(2)) * SDSD #measures the width of poincare cloud
    SD2 = np.sqrt((2 * SDRR ** 2) - (0.5 * SDSD ** 2)) #measures the length of the poincare cloud
 
    return {'SD1': SD1, 'SD2': SD2}
    
def hraMethod(RRints):
    """
    Perform analysis to quantify heart rate assymmetry (HRA).
    
    Input           :
    
     - RRints       : [list] of RR intervals
        
    Output          :

     - C_DOWN, C_UP : {dict} with keys 'C_DOWN' (numpy.float64) and 
                      'C_UP' (numpy.float64)    
     """    
    ax1 = RRints[:-1]
    ax2 = RRints[1:]
    SD1I = np.sqrt((1 / len(ax1)) * (np.sum((ax1 - ax2) ** 2) / 2))
    ax1ax2 = (ax1 - ax2) / np.sqrt(2)
    indices_up = np.where(ax1ax2 > 0)
    indices_down = np.where(ax1ax2 < 0)
    SD1_UP = np.sqrt(np.sum(ax1ax2[indices_up] ** 2) / len(ax1))
    SD1_DOWN = np.sqrt(np.sum(ax1ax2[indices_down] ** 2) / len(ax1))    
    C_UP = SD1_UP ** 2 / SD1I ** 2 #for decelerations
    C_DOWN = SD1_DOWN ** 2 / SD1I ** 2 #for accelerations
    
    return {'C_UP': C_UP, 'C_DOWN': C_DOWN}

def correlationCoef(RRints):
    """
    Computes interbeat autocorrelation coefficient
    
    Input    :
    
     - RRints: [list] of RR intervals
        
    Output   : 
        
     - r_rr  : [numpy.float64] interbeat autocorrelation coefficient
 
     """    
    ax1 = RRints[:-1]
    ax2 = RRints[1:]
    mu_rr = np.mean(RRints)
    r_rr = np.mean((ax1 - mu_rr) * (ax2 - mu_rr)) / (np.sqrt(np.mean((ax1 - mu_rr) ** 2) * np.mean((ax2 - mu_rr) ** 2)))

    return r_rr
    
def plotRRintHist(RRints):
    """ 
    Histogram distribution of poincare points projected onto the x-axis
    
    Input    :
    
     - RRints: [list] of RR intervals
        
    Output   :
        
     - RR interval histogram plot    
    """    
    plt.hist(RRints, bins = 'auto')
    plt.xlabel('RR Interval')
    plt.ylabel('Number of RR Intervals')
    plt.title('RR Interval Histogram')
    plt.show()

def plotWidthHist(RRints):    
    """  
    Histogram distribution of poincare points projected along the direction of 
    line-of-identity, or along the line perpendicular to the line-of-identity.
    
    Input    :
    
     - RRints: [list] of RR intervals
        
    Output   :
        
     - 'Width', or delta-RR interval, histogram plot      
     """   
    ax1 = RRints[:-1]
    ax2 = RRints[1:]
    x1 = (np.cos(np.pi / 4) * ax1) - (np.sin(np.pi / 4) * ax2)
    plt.hist(x1, bins = 'auto')
    plt.title('Width (Delta-RR Interval) Histogram')
    plt.show()
    
def plotLengthHist(RRints):    
    """
    Histogram distribution of poincare points projected along the line-of-identty.
    
    Input    :
    
     - RRints: [list] of RR intervals
        
    Output   :
        
     - 'Length' histogram plot
     """
     
    ax1 = RRints[:-1]
    ax2 = RRints[1:]
    x2 = (np.sin(np.pi / 4) * ax1) + (np.cos(np.pi / 4) * ax2)
    plt.hist(x2, bins = 'auto')
    plt.title('Length Histogram')
    plt.show()






