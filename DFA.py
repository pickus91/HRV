# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 14:35:26 2017

@author: picku
"""
import numpy as np
import matplotlib.pyplot as plt

      
def dfa(timeSeries, scale, q, m):
    #m: [int] order of polynomial fit
    #q: [int] order of fluctuation coefficient        
    
    #Step(1): Determine the "profile" (integrated signal with subtracted offset)
    integratedSeries = np.cumsum(timeSeries - np.mean(timeSeries)) #y_k
    
    #Step(2): Divide profile into N non-overlapping segments of equal length s
    shape = (integratedSeries.shape[0] // scale, scale)
    nwLength = shape[0] * shape[1] 
    windowedData = np.reshape(integratedSeries[0:nwLength], shape)    
    
    #Step(3): Calculate local trend for each 2Ns segments by a least squares fit of the series. Then determine the variance for each segment v, v = 1,...,Ns
    x = np.arange(windowedData.shape[1]) 
    rms = np.empty(windowedData.shape[0])  
    
    for i, segment in enumerate(windowedData):
        
        pl = np.polyfit(x, segment, m) 
        y_n = np.polyval(pl, x)
        rms[i] = np.sqrt(np.mean((segment - y_n) ** 2)) 
      
    #(Step 4): Average over all segments to obtain the qth order fluctuation coefficient      
    F = np.mean(rms ** q) ** (1 / q)
    
    return F

def scalingExponent(timeSeries, lowerScaleLimit, upperScaleLimit, scaleDensity = 100, m, q, plot = 0):
    
    startBeatNum = np.log(startBeatNum) / np.log(10)
    stopBeatNum = np.log(stopBeatNum) / np.log(10)
    
    scales = np.floor(np.logspace(np.log10(10 ** startBeatNum), np.log10(10 ** stopBeatNum), scaleDensity))
    
    F = np.zeros(scales.shape[0])  
    
    for j,scale in enumerate(scales):
        F[j] = dfa(timeSeries, scale, q, m)  #timeSeries = RR series to be analyzed
        
    pl2 = np.polyfit(np.log2(scales), np.log2(F), 1) #linear fit   
    lineFit = np.polyval(pl2, np.log2(scales))
    
    scaleParameter = pl2[0] #Finding scaling exponent (Hurst exponent = h(m = 2))
    
    if plot == 1:
        
        plt.loglog(scales, F, '.', color = c , markersize = 3)
        plt.plot(scales, 2 ** lineFit, color = c, linewidth = 1, label = '%s = %0.2f' % (l,scaleParameter))
        plt.xlabel('log10(scale)')
        plt.ylabel('log10(F)')
        plt.title('F vs Scale')
        plt.legend(loc='best')
        
        
        
        
    

        
        





