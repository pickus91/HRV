# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 14:03:46 2017

@author: picku

multiScaleEntropy implements the algorithm to calculate multiscale entropy (MSE) for
a complex time series as shown in:
    
Costa, Madalena, Ary L. Goldberger, and C-K. Peng. "Multiscale entropy analysis of 
complex physiologic time series." Physical review letters 89.6 (2002): 068102
                                
The MSE algorithm works by computing sample entropy (SampEn) at multiple scales, making
it advantageous for analyzing features related to structure on scales other than the
shortest one. Sample entropy is a modification of approximate entropy (ApEn) in that 
(1) it does not count self-matching and (2) it does not depend as much on the 
length of the time series. Sample entropy algorithm details can be found in:
    
Pincus, Steven M. "Approximate entropy as a measure of system complexity." 
Proceedings of the National Academy of Sciences, 88(6) (1991): 2297-2301
"""
import numpy as np

#For sleep staging analysis
n = 30
startBeatNum = 10
sBn = np.log(startBeatNum) / np.log(10)
stopBeatNum = 40
stBn = np.log(stopBeatNum) / np.log(10)
scales = np.floor(np.logspace(np.log10(10 ** sBn), np.log10(10 ** stBn), n))

def multiScaleEntropy(timeSeries, scales, r, m):
    '''MultiScaleEntropy - algorithm to calculate multiscale entropy (MSE) for complex time series.
       Introduced by Costa, et al (2002)
    
    Inputs Arguments:         
        - timeSeries: [list] of time series data 
        - scales    : [list] of scales for MSE calculation
        - r         : [float] tolerance (typically 0.2*std(time series))
        - m         : [int] embedding dimension (typically 2)

    Ouputs Arguments:
        - MSE       : [array] of sample entropy at scales --> size [1 length(scales)]
    '''
    N = len(timeSeries)
    MSE = []
    
    for s in scales:        
        if s > N:
            print('Scale must not exceed length of the time series')
        else:
            #"course-graining" process      
            cuts = np.reshape(timeSeries[:(N - int(np.mod(N, s)))], (int(np.floor(N/s)), int(s)))
            coarseGrainedSeries = np.mean(cuts, 1)
            MSE.append(sampEn(coarseGrainedSeries, m, r))       
            
    return MSE        
        
def sampEn(x, m, r): 
    '''sampEn computes sample entropy of time series 'x'. 
        
        Input Arguments:
            - x   : [list] of time-series data 
            - m   : [int] embedding dimension (typically 2)
            - r   : [float] tolerance (typically 0.2*std(x))
            
        Output Arguments:
            - sE  : [float] sample entropy value
    '''
    correlHolder = []
    L = len(x)
    W = [m, m + 1]
   
    for i in W: 
        #create space for extracted windows
        z = np.zeros((L - i + 1 , i))
        #extract windows and store in matrix z       
        j = np.arange(0, L - i + 1)       
        for k in j:
            z[k][:] = x[k:k + i]      
       
        count = 0 
        for l in j:
            for n in j:
                if n != l:
                   #Chebyshev distance
                    D = max(abs(z[l][:] - z[n][:]))
                    if D < r: #distance exceeds tolerance
                        count += 1                       
        
        correlHolder.append(count)
        
    B = correlHolder[0]  
    A = correlHolder[1]
              
    if (A == 0) | (B == 0): #need to do some research on this result***
        return 0
    else:    
        return -1 * np.log(A / B)
   
   
   
   

