# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 12:08:40 2017

@author: picku
"""
import numpy as np
from scipy import interpolate, signal

def frequencyDomain(RRints):

    #Remove ectopic beats
    NNs = removeBadRRs(RRints)
    
    #Resample @ 4 Hz
    fsResamp = 4   
    tmStamps = np.cumsum(NNs) #in seconds 
    f = interpolate.interp1d(tmStamps, NNs, 'cubic')
    tmInterp = np.arange(tmStamps[0], tmStamps[-1], 1/fsResamp)
    RRinterp = f(tmInterp)          
    
    #Remove DC component     
    RRseries = RRinterp - np.mean(RRinterp)
        
    #Pwelch w/ zero pad for frequecy resolution of approx 0.0001    
    fxx, pxx = signal.welch(RRseries, fsResamp, nfft = 1024, window = 'hann')
    '''
    plt.plot(fxx,pxx,color = 'b')
    plt.xlim([0,0.5])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('PSD(ms^2/Hz)')
    plt.title('RR Spectrum')
    '''
    df = fxx[1] - fxx[0]
    vlf= (0.003, 0.04)
    lf = (0.04, 0.15)
    hf = (0.15, 0.4)

    vlf_power = np.trapz(pxx[np.logical_and(fxx >= vlf[0], fxx < vlf[1])], dx = df)      
    lf_power = np.trapz(pxx[np.logical_and(fxx >= lf[0], fxx < lf[1])], dx = df)            
    hf_power = np.trapz(pxx[np.logical_and(fxx >= hf[0], fxx < hf[1])], dx = df)             
    totalPower = vlf_power + lf_power + hf_power
    
    #Normalize and take log
    vlf_NU_log = np.log((vlf_power / (totalPower - vlf_power)) + 1)
    lf_NU_log = np.log((lf_power / (totalPower - vlf_power)) + 1)
    hf_NU_log = np.log((hf_power / (totalPower - vlf_power)) + 1)
    lfhfRation_log = np.log((lf_power / hf_power) + 1)   
    
    freqDomainFeats = {'VLF_Power': vlf_NU_log, 'LF_Power': lf_NU_log,
                       'HF_Power': hf_NU_log, 'LF/HF': lfhfRation_log}
    
    return freqDomainFeats
        
def removeBadRRs(RRints):
    #RR intervals differing by more than 20% from the one proceeding it are removed
    NNs = []
    for c, rr in enumerate(RRints):        
        if abs(rr - RRints[c-1]) <= 0.20 * RRints[c-1]:
            NNs.append(rr)
            
    return NNs
            
    
    
    
    
          

    
    
    