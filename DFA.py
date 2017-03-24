# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 14:35:26 2017

@author: picku

This code carries out multifractal detrended fluctuation analysis (MF-DFA) as
described in:

Kantelhardt, Jan W., et al. "Multifractal detrended fluctuation analysis of 
nonstationary timer series." Physica A: Statistical Mechanics and its Applications
316.1 (2002): 87-114
    
Unlike standard DFA, which is used to analyze (mono-) fractal scaling properties 
of nonstationary time series, MF-DFA allows for analysis of multi-fractal
nonstationary time series using one additional step - a q dependent averaging
procedure - where q = 2 performs the standard DFA procedure developed in:
    
Peng, C-K., et al. "Quantification of scaling exponents and crossover phenomena
in nonstationary heartbeat time series." Chaos: An Interdisciplinary Journal of 
Nonlinear Science 5.1 (1995): 82-87
    
Characterizing the short- and long-range correrlation properties of the RR interval
time series (a non-stationary physiological signal) is helpful in determining
fluctuations due to the intrinsic dynamics of the cardio-respiratory system itself
and those fluctuations that arise from external (environmental) stimuli. The scaling 
exponent, alpha, can be used to represent the autocorrelatin properties of the signal:

alpha < 1/2 : anti-correlated
alpha = 1/2 : uncorrelated (white noise)
alpha > 1/2 : correlated
alpha = 1   : 1/f-noise
alpha > 1   : non-stationary
alpha = 3/2 : Brownian noise (random walk)
                                      
"""
import numpy as np
import matplotlib.pyplot as plt
      
def dfa(timeSeries, scale, q, m):
    '''
    Input Arguments  :
        
        -timeSeries  : [list] of RR intervals
        -scale       : [int]
        -m           : [int] order of polynomial fit
        -q           : [int] order of fluctuation coefficient
        
    Output Arguments:
        
        -F           : [float] fluctuation coefficient
    
    '''  
    if len(timeSeries) < scale:
        print('Length of time series must be greater than scale')
    
    else:    
        #Step(1): Determine the "profile" (integrated signal with subtracted offset)
        integratedSeries = np.cumsum(timeSeries - np.mean(timeSeries)) #y_k
        
        #Step(2): Divide profile into N non-overlapping segments of equal length s
        shape = (int(integratedSeries.shape[0]) // int(scale), int(scale))
        nwLength = shape[0] * shape[1] 
        windowedData = np.reshape(integratedSeries[0:int(nwLength)], shape)    
        
        #repeate same procedure from opposite end of integrated series
        windowedData2 = np.reshape(integratedSeries[len(integratedSeries) - int(nwLength):], shape)

        segments = np.concatenate((windowedData, windowedData2))
        #Step(3): Calculate local trend for each 2Ns segments by a least squares fit of the series. 
        #Then determine the variance for each segment v, v = 1,...,Ns
        x = np.arange(segments.shape[1]) 
        rms = np.empty(segments.shape[0])  
        
        for i, segment in enumerate(segments):
            
            pl = np.polyfit(x, segment, m) 
            y_n = np.polyval(pl, x)
            rms[i] = np.sqrt(np.mean((segment - y_n) ** 2)) 
          
        #(Step 4): Average over all segments to obtain the qth order fluctuation coefficient      
        F = np.mean(rms ** q) ** (1 / q)
    
        return F

def scalingExponent(timeSeries, lowerScaleLimit, upperScaleLimit, scaleDensity, m, q, plot):
    '''
    Input Arguments   :
        
    - timeSeries      : [list] of RR intervals
    - lowerScaleLimit : [int] minimum of scales for MF-DFA
    - upperScaleLimit : [int] maximum of scales for MF-DFA
    - scaleDensity    : [int] number of scales between lowerScaleLimit and upperScaleLimit in which DFA is conducted
    - m               : [int] order of polynomial fit
    - q               : [int] order of fluctuation coefficient
    - plot            : [1] --> plot of log(scales) vs log(F)
        
    Output Arguments  :
        
    - H               : [float] scaling exponent 
    
    '''
    startBeatNum = np.log(lowerScaleLimit) / np.log(10)
    stopBeatNum = np.log(upperScaleLimit) / np.log(10)
    
    scales = np.floor(np.logspace(np.log10(10 ** startBeatNum), np.log10(10 ** stopBeatNum), scaleDensity))
    
    F = np.zeros(scales.shape[0])  
    
    for j,scale in enumerate(scales):
        F[j] = dfa(timeSeries, int(scale), q, m)  #timeSeries = RR series to be analyzed
      
    #Step(5) Determine scale behavior of the fluctuation functions by analyzing
    #log-log plots of F versus s for each value of q
    pl2 = np.polyfit(np.log2(scales), np.log2(F), 1) #m = 1 ---> linear fit   
    lineFit = np.polyval(pl2, np.log2(scales))
    
    scaleParameter = pl2[0] #Finding scaling exponent (Hurst exponent = h(m = 2))
    
    if plot == 1:       
        
        plt.loglog(scales, F, '.', markersize = 3)
        plt.plot(scales, 2 ** lineFit, linewidth = 1, label = r'$\alpha$ = %0.2f' % (scaleParameter))
        plt.xlabel(r'$\log_{10}$(scale)')
        plt.ylabel(r'$\log_{10}$(F)')
        plt.title('F vs Scale')
        plt.legend(loc='best')
        
        



    

        
        





