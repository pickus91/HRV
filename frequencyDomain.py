# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 12:08:40 2017

@author: picku

"""
import numpy as np
from scipy import interpolate, signal
import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')
import matplotlib.patches as mpatches
from collections import OrderedDict

def frequencyDomain(RRints, band_type = None, lf_bw = 0.11, hf_bw = 0.1, plot = 0):
    """ Computes frequency domain features on RR interval data
    
    Parameters:
    ------------
    RRints : list, shape = [n_samples,]
           RR interval data
    
    band_type : string, optional
             If band_type = None, the traditional frequency bands are used to compute 
             spectral power:
           
                 LF: 0.003 - 0.04 Hz
                 HF: 0.04 - 0.15 Hz
                 VLF: 0.15 - 0.4 Hz           
           
             If band_type is set to 'adapted', the bands are adjusted according to 
             the protocol laid out in:
           
             Long, Xi, et al. "Spectral boundary adaptation on heart rate 
             variability for sleep and wake classification." International 
             Journal on Artificial Intelligence Tools 23.03 (2014): 1460002. 
                        
    lf_bw : float, optional
          Low frequency bandwidth centered around LF band peak frequency
          when band_type is set to 'adapted'. Defaults to 0.11
             
    hf_bw : float, optional
          High frequency bandwidth centered around HF band peak frequency
          when band_type is set to 'adapted'. Defaults to 0.1
          
    plot : int, 1|0
          Setting plot to 1 creates a matplotlib figure showing frequency
          versus spectral power with color shading to indicate the VLF, LF,
          and HF band bounds.
    
    Returns:
    ---------
    freqDomainFeats : dict
                   VLF_Power, LF_Power, HF_Power, LF/HF Ratio              
    """

    #Remove ectopic beats
    #RR intervals differing by more than 20% from the one proceeding it are removed
    NNs = []
    for c, rr in enumerate(RRints):        
        if abs(rr - RRints[c-1]) <= 0.20 * RRints[c-1]:
            NNs.append(rr)
      
    #Resample @ 4 Hz
    fsResamp = 4   
    tmStamps = np.cumsum(NNs) #in seconds 
    f = interpolate.interp1d(tmStamps, NNs, 'cubic')
    tmInterp = np.arange(tmStamps[0], tmStamps[-1], 1/fsResamp)
    RRinterp = f(tmInterp)          
    
    #Remove DC component     
    RRseries = RRinterp - np.mean(RRinterp)
        
    #Pwelch w/ zero pad     
    fxx, pxx = signal.welch(RRseries, fsResamp, nfft = 2**14, window = 'hann')    
    
    vlf= (0.003, 0.04)
    lf = (0.04, 0.15)
    hf = (0.15, 0.4)
    
    plot_labels = ['VLF', 'LF', 'HF']
        
    if band_type == 'adapted':     
            
        vlf_peak = fxx[np.where(pxx == np.max(pxx[np.logical_and(fxx >= vlf[0], fxx < vlf[1])]))[0][0]] 
        lf_peak = fxx[np.where(pxx == np.max(pxx[np.logical_and(fxx >= lf[0], fxx < lf[1])]))[0][0]]
        hf_peak = fxx[np.where(pxx == np.max(pxx[np.logical_and(fxx >= hf[0], fxx < hf[1])]))[0][0]]
    
        peak_freqs =  (vlf_peak, lf_peak, hf_peak) 
            
        hf = (peak_freqs[2] - hf_bw/2, peak_freqs[2] + hf_bw/2)
        lf = (peak_freqs[1] - lf_bw/2, peak_freqs[1] + lf_bw/2)   
        vlf = (0.003, lf[0])
        
        if lf[0] < 0:
            print('***Warning***: Adapted LF band lower bound spills into negative frequency range')
            print('Lower thresold of LF band has been set to zero')
            print('Adjust LF and HF bandwidths accordingly')
            lf = (0, lf[1])        
            vlf = (0, 0)
        elif hf[0] < 0:
            print('***Warning***: Adapted HF band lower bound spills into negative frequency range')
            print('Lower thresold of HF band has been set to zero')
            print('Adjust LF and HF bandwidths accordingly')
            hf = (0, hf[1])        
            lf = (0, 0)        
            vlf = (0, 0)
            
        plot_labels = ['Adapted_VLF', 'Adapted_LF', 'Adapted_HF']

    df = fxx[1] - fxx[0]
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
                       
    if plot == 1:
        #Plot option
        freq_bands = {'vlf': vlf, 'lf': lf, 'hf': hf}
        freq_bands = OrderedDict(sorted(freq_bands.items(), key=lambda t: t[0]))
        colors = ['lightsalmon', 'lightsteelblue', 'darkseagreen']
        fig, ax = plt.subplots(1)
        ax.plot(fxx, pxx, c = 'grey')
        plt.xlim([0, 0.40])
        plt.xlabel(r'Frequency $(Hz)$')
        plt.ylabel(r'PSD $(s^2/Hz$)')
        
        for c, key in enumerate(freq_bands):
            ax.fill_between(fxx[min(np.where(fxx >= freq_bands[key][0])[0]): max(np.where(fxx <= freq_bands[key][1])[0])],
                            pxx[min(np.where(fxx >= freq_bands[key][0])[0]): max(np.where(fxx <= freq_bands[key][1])[0])],
                            0, facecolor = colors[c])
            
        patch1 = mpatches.Patch(color = colors[0], label = plot_labels[2])
        patch2 = mpatches.Patch(color = colors[1], label = plot_labels[1])
        patch3 = mpatches.Patch(color = colors[2], label = plot_labels[0])
        plt.legend(handles = [patch1, patch2, patch3])
        plt.show()

    return freqDomainFeats
        
            
    
    
    
    
          

    
    
    
