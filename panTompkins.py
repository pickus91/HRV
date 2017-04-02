# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 17:00:49 2017

@author: picku

This function implements the Pan-Tomkins ECG QRS detection algorithm
as described in:
    
Pan, Jiapu, and Willis J. Tompkins. "A real-time QRS detection algorithm." 
IEEE transactions on biomedical engineering 3 (1985): 230-236.

The algorithm utilizes filtering, adaptive thresholding, and criteria based on
human cardiac physiology to detect QRS complexes in the face of noise and quickly 
changing and diverse ECG morphologies. This function implements the Pan-Tompkins
algorithm as it was originally published, along with two modifications which
include additional filtering and eliminating potential QRS detections that occur
within the refractory period. Since this algorithm is often used to find R-peak
locations (and not just general QRS detection) for applications such as Heart
Rate Variability (HRV) analysis, this function also performs a neighborhood search
around the final QRS detection locations to find exact R-peak locations. 
A summary of the algorithm is described below:

(1) Band Pass Filtering

    The ECG signal is first bandpassed filtered between 5 and 15 Hz to eliminate
    noise due to muscle contractions, 60 Hz noise, baseline wander, and T-wave
    interference. As implemented in the Pan-Tompkins publication, a low and 
    high pass filter are applied in series to acheive a frequency passband between
    5 and 12 Hz. 
    
    LP ---> y(nT) = 2y(nT-T) - y(nT-2T) + x(nT) - 2x(nT - 6) + x(nT-12T)
    
    HP ---> y(nT) = 32x(nT - 16T) - [y(nT - T) + x(nT) - x(nT - 32T)]
    
    This code also applies an additional 5th order Butterworth filter with 
    cutoff frequencies of 5 and 15 Hz to make this algorithm robust to noisier 
    data sets.                     

(2) Derivative

    QRS slope information is computed via differentiation of the filtered ECG
    signal.
    
           y(nT) = (1/8T)[-x(nT-2T) - 2x(nT - T) + 2x(nT +T) + x(nT + 2T)]
                              
(3) Squaring Function

    The squaring function serves for non-linear amplification the high frequency 
    components of the derivative signal associated with the QRS complexes.
    
                             y(nT) = [x(nT)] ** 2
           
(4) Moving Integration Waveform

    A final waveform is created using a moving-window integrator. This waveform
    serves to obtain waveform feature information, with each rising edge
    corresponding to the location of the QRS complex (increased amplitude --> increased 
    area)
    
           y(nT) = (1/N)[x(nT - (N - 1)T + x(nT - (N - 2)T) + ... + x(nT)]
    
    The width of the waveform was found empiracally in the publication to be
    30 samples wide for a sample frequency of 200 Hz. In order to increase
    the flexibility of this algorithm to signals sampled at different 
    frequencies, the window width is automatically adjusted so that the window
    is approximately 150 ms wide regardless of sampling frequency.     
    
(5) Fudicial Mark

    The rising edge of the integration waveform corresponds with the QRS complex.
    A fudicial mark indicating the temporal location of the QRS complex is made by 
    locating the peaks of the integration waveform associated with the R peak of 
    the QRS complex. This is done through differentiation, zero-crossing location, 
    and moving average amplitude thresholding to eliminate peaks from the P and T
    wave features.  
    
(6) Adjusting Signal and Noise Thresholds

     Following a brief two second initilization phase, two sets of signal and noise 
     thresholds are adjustedfor both the integration waveform and the filtered ECG
     signal, respectively. These thresholds are based off previous peaks determined 
     to be either signal peaks or noise peaks. Thresholds are able to adapt quickly to 
     rapidly changing heart rates by keeping running estimates of both signal and 
     noise levels from previous peak assignments. Thresholds are computed as follows
     for both the integration waveform and the filtered ECG:
     
     If signal peak  ---> SPK = 0.125 * Peak_Amplitude + 0.875 * SPK
     If noise peak   ---> NPK = 0.125 * Peak_Amplitude + 0.875 * NKP
     
                          THRESHOLD1 = NPK + 0.25 * (SPK - NPK)
                          THRESHOLD2 = 0.5 * THRESHOLD1
        
     A peak is considered to be a signal peak if it exceeds THRESHOLD1 or
     exceeds THRESHOLD2 if a searchback is triggered. 
        
(6) RR Rate Limits
    
    Two heartbeats are needed to establish the average RR interval and rate
    limits. If any of the eight most recent sequential beats fall outside the the 
    accepted low and high RR-interval limits, heart rate is considered to be 
    irregular and the signal and noise thresholds are reduced by half in order to
    increase sensitivity.
    
(6) Searchback

    If a QRS complex is not found within 166% of the average of the previous
    eight beats during normal sinus rhythm, it is assumed that a QRS complex
    has been missed. A searchback process is then triggered, which finds 
    the maximal peak within the current signal and noise thresholds to be
    a QRS complex candidate. If this QRS complex candidate exceeds THRESHOLD2,
    the signal levels will take more consideration of the current peak 
    amplitude and less on the previous signal values via the following modification
    to the signal level running estimate:
        
                    SPK = 0.25 * Peak_Amplitude + 0.75 * SPKI

(7) T-wave Identification

    Once a QRS complex is identified, there is a 200 ms refractory in which it is 
    physiologically impossible for another beat to occur, thus allowing for the 
    elimination of any QRS detection during this time frame. If a QRS is detected
    after the refractory period but before 360 ms after the temporal location
    of the previous QRS detection, we must decide whether this peak is an actual
    QRS complex or a T-wave. If the slope of this peak is less than half of the 
    slope of the previous QRS, it is identified as a T-wave.

(8) Detection

    A QRS is identified if it is identified in both the band-passed filtered ECG signal 
    and in the integration waveform. 
    
(9) Neighborhood Search

    A neighborhood search with a search window of 0.2 seconds is used to find the 
    highest amplitude point in the general vicinity of the QRS complex. This is 
    the location of the R-peak.

"""
import numpy as np
from matplotlib import style
from scipy import signal
import matplotlib.pyplot as plt 
style.use('ggplot')

def panTompkins(ECG, fs, plot = 1): 
    """
    Inputs:
    - ECG   : [list] | [numpy.ndarray] of ECG samples
    - fs    : [int] sampling frequency
    - plot  : [1|0] optional plot of R-peak detections overlayed on ECG signal

    Outputs:
    - Rpeaks : [list] of integers indicating R peak sample number locations
    """    
    if type(ECG) == list or type(ECG) is np.ndarray:
        ECG = np.array(ECG)             
        
    #Initialize
    RRAVERAGE1 = []
    RRAVERAGE2 = []
    IWF_signal_peaks = []
    IWF_noise_peaks = []
    noise_peaks = []
    ECG_bp_peaks = np.array([])
    ECG_bp_signal_peaks = []
    ECG_bp_noise_peaks = []
    final_R_locs = []
    T_wave_found = 0      
    
    #LOW PASS FILTERING
    #Transfer function: H(z)=(1-z^-6)^2/(1-z^-1)^2
    a = np.array([1, -2, 1])
    b = np.array([1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1])   
        
    impulse = np.repeat(0., len(b)); impulse[0] = 1.    
    impulse_response = signal.lfilter(b,a,impulse)
    
    #convolve ECG signal with impulse response
    ECG_lp = np.convolve(impulse_response, ECG)
    ECG_lp = ECG_lp / (max(abs(ECG_lp)))
    delay = 12 #full convolution
    
    #HIGH PASS FILTERING
    #Transfer function: H(z)=(-1+32z^-16+z^-32)/(1+z^-1)
    a = np.array([1, -1])           
    b = np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 32, -32, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, -1])
                  
    impulse = np.repeat(0., len(b)); impulse[0] = 1.    
    impulse_response = signal.lfilter(b,a,impulse)
    
    ECG_lp_hp = np.convolve(impulse_response, ECG_lp)
    ECG_lp_hp = ECG_lp_hp/(max(abs(ECG_lp_hp)))
    delay = delay + 32 
    
    #BAND PASS FILTER 
    nyq = fs / 2        
    lowCut = 5 / nyq  #cut off frequencies are normalized from 0 to 1, where 1 is the Nyquist frequency
    highCut = 15 / nyq
    order = 5
    b,a = signal.butter(order, [lowCut, highCut], btype = 'bandpass')
    ECG_bp = signal.lfilter(b, a, ECG_lp_hp)
    
    #DIFFERENTIATION
    #Transfer function: H(z)=(1/8T)(-z^-2-2z^-1+2z^1+z^2)
    T = 1/fs
    b = np.array([-1, -2, 0, 2, 1]) * (1 / (8 * T))
    a = 1
    #Note impulse response of the filter with a = [1] is b
    ECG_deriv = np.convolve(ECG_bp, b)
    delay = delay + 4 
    
    #SQUARING FUNCTION
    ECG_squared = ECG_deriv ** 2
    
    #MOVING INTEGRATION WAVEFORM 
    N = int(np.ceil(0.150 * fs)) 
    ECG_movavg = np.convolve(ECG_squared,(1 / N) * np.ones((1, N))[0])
    
    #FUDICIAL MARK ON MOVING INTEGRATION WAVEFORM
    peaks = findPeaks(ECG_movavg)
    
    #LEARNING PHASE 1   
    #2 second initialize phase for MIW, 25% of max amplitude considered signal, 50% of mean signal considered noise
    initializeTime = 2 * fs 
    SPKI = max(ECG_movavg[:initializeTime]) * 0.25 
    NPKI = np.mean(ECG_movavg[:initializeTime]) * 0.5 
    THRESHOLDI1 = NPKI + 0.25 * (SPKI-NPKI)
    THRESHOLDI2 = 0.5 * THRESHOLDI1 
    
    #2 second initialize for filtered signal, 25% of max amplitude considered signal, 50% of mean signal considered noise
    initializeTime = 2 * fs 
    SPKF = max(ECG_bp[:initializeTime]) * 0.25 
    NPKF = np.mean(ECG_bp[:initializeTime]) * 0.5 
    THRESHOLDF1 = NPKF + 0.25 * (SPKF-NPKF)
    THRESHOLDF2 = 0.5 * THRESHOLDF1
    
    peaks = peaks[peaks > initializeTime] #ignore peaks that occur during initialization window

    for c,peak in enumerate(peaks):                
        #find corresponding peaks in filtered ECG using neighborhood search window +- 0.15 seconds       
        searchInterval = int(np.round(0.15 * fs))
        searchIndices = np.arange(peak - searchInterval, peak + searchInterval + 1, 1)        
        #neighborhood search indices cannot be negative and cannot exceed length of filtered ECG
        if searchIndices[0] >= 0 and all(searchIndices <= len(ECG_bp)):            
             ECG_bp_peaks = np.append(ECG_bp_peaks, np.where(ECG_bp == max(ECG_bp[searchIndices]))[0][0])          
        else:
             ECG_bp_peaks = np.append(ECG_bp_peaks, np.where(ECG_bp == max(ECG_bp[searchIndices[0]:len(ECG_bp)-1])))
        #LEARNING PHASE 2
        if c > 0 and c < len(ECG_bp_peaks):                     
            if c < 8:                
                RRAVERAGE1_vec = np.diff(peaks[:c + 1]) / fs
                RRAVERAGE1_mean = np.mean(RRAVERAGE1_vec)
                RRAVERAGE1.append(RRAVERAGE1_mean) 
                
                RR_LOW_LIMIT = 0.92 * RRAVERAGE1_mean
                RR_HIGH_LIMIT = 1.16 * RRAVERAGE1_mean
                RR_MISSED_LIMIT = 1.66 * RRAVERAGE1_mean                   
            else:                
                RRAVERAGE1_vec = np.diff(peaks[c - 8:c + 1]) / fs
                RRAVERAGE1_mean = np.mean(RRAVERAGE1_vec)
                RRAVERAGE1.append(RRAVERAGE1_mean) 
    
                for rr in np.arange(0, len(RRAVERAGE1_vec)):
                    if RRAVERAGE1_vec[rr] > RR_LOW_LIMIT and RRAVERAGE1_vec[rr] < RR_HIGH_LIMIT:                              
                        RRAVERAGE2.append(RRAVERAGE1_vec[rr])                                     
                        if len(RRAVERAGE2) > 8:
                            del RRAVERAGE2[:len(RRAVERAGE2) - 8]
    
                if len(RRAVERAGE2) == 8:
                    RR_LOW_LIMIT = 0.92 * np.mean(RRAVERAGE2)        
                    RR_HIGH_LIMIT = 1.16 * np.mean(RRAVERAGE2)
                    RR_MISSED_LIMIT = 1.66 * np.mean(RRAVERAGE2)
            #If irregular heart beat detected in previous 9 beats, lower signal thresholds by half to increase detection sensitivity            
            current_RR_movavg = RRAVERAGE1[-1] 
            if current_RR_movavg < RR_LOW_LIMIT or current_RR_movavg > RR_MISSED_LIMIT: 
                #MIW thresholds        
                THRESHOLDI1 = 0.5 * THRESHOLDI1
                THRESHOLDI2 = 0.5 * THRESHOLDI1
                #Filtered ECG thresholds
                THRESHOLDF1 = 0.5 * THRESHOLDF1
                THRESHOLDF2 = 0.5 * THRESHOLDF1
               
            #Search back triggered if current RR interval is greater than RR_MISSED_LIMIT
            currentRRint = RRAVERAGE1_vec[-1]
            if currentRRint > RR_MISSED_LIMIT:  
                SBinterval = int(np.round(currentRRint * fs))
                #find local maximum in the search back interval between signal and noise thresholds                        
                SBdata_IWF = ECG_movavg[peak - SBinterval + 1:peak + 1]               
                
                SBdata_IWF_filtered = np.where((SBdata_IWF > THRESHOLDI1))[0]
                SBdata_max_loc = np.where(SBdata_IWF == max(SBdata_IWF[SBdata_IWF_filtered]))[0][0]

                if len(SBdata_IWF_filtered) > 0:   
                    SB_IWF_loc = peak - SBinterval + 1 + SBdata_max_loc
                    IWF_signal_peaks.append(SB_IWF_loc) 
                    #update signal and noise thresholds
                    SPKI = 0.25 * ECG_movavg[SB_IWF_loc] + 0.75 * SPKI                         
                    THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI)
                    THRESHOLDI2 = 0.5 * THRESHOLDI1               
                    #finding corresponding search back peak in ECG bandpass using 0.15 s neighborhood search window
                    if SB_IWF_loc < len(ECG_bp):
                        SBdata_ECGfilt = ECG_bp[SB_IWF_loc - round(0.15 * fs): SB_IWF_loc]                    
                        SBdata_ECGfilt_filtered = np.where((SBdata_ECGfilt > THRESHOLDF1))[0]
                        SBdata_max_loc2 = np.where(SBdata_ECGfilt == max(SBdata_ECGfilt[SBdata_ECGfilt_filtered]))[0][0]
                                     
                    else:
                        SBdata_ECGfilt = ECG_bp[SB_IWF_loc - round(0.15 * fs):]
                        SBdata_ECGfilt_filtered = np.where((SBdata_ECGfilt > THRESHOLDF1))[0]
                        SBdata_max_loc2 = np.where(SBdata_ECGfilt == max(SBdata_ECGfilt[SBdata_ECGfilt_filtered]))[0][0]

                            
                    if ECG_bp[SB_IWF_loc - round(0.15 * fs) + SBdata_max_loc2] > THRESHOLDF2: #QRS complex detected in filtered ECG
                        #update signal and noise thresholds                                                          
                        SPKF = 0.25 * ECG_bp[SB_IWF_loc - round(0.15 * fs) + SBdata_max_loc2] + 0.75 * SPKF                            
                        THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF)
                        THRESHOLDF2 = 0.5 * THRESHOLDF1                            
                        ECG_bp_signal_peaks.append(SB_IWF_loc - round(0.15 * fs) + SBdata_max_loc2)                                                 
    
            #T-WAVE AND QRS DISRCIMINATION    
            if ECG_movavg[peak] >= THRESHOLDI1: 
                if currentRRint > 0.20 and currentRRint < 0.36 and c > 0: 
                    #Slope of current waveform (possible T wave)
                    #mean width of QRS complex: 0.06 - 0.10 sec         
                    maxSlope_current = max(np.diff(ECG_movavg[peak - round(fs * 0.075):peak + 1]))
                    #slope of the waveform (most likely QRS) that preceeded it
                    maxSlope_past = max(np.diff(ECG_movavg[peaks[c - 1] - round(fs * 0.075): peaks[c - 1] + 1]))
                    if maxSlope_current < 0.5 * maxSlope_past: #T-wave found                        
                        T_wave_found = 1                
                        #keep track of peaks marked as 'noise'
                        IWF_noise_peaks.append(peak)                
                        #update Noise levels
                        NPKI = 0.125 * ECG_movavg[peak] + 0.875 * NPKI                                            
                               
                if not T_wave_found: #current peak is a signal peak                    
                    IWF_signal_peaks.append(peak)
                    #adjust signal levels
                    SPKI = 0.125 * ECG_movavg[peak]  + 0.875 * SPKI
                    #check if corresponding peak in filtered ECG is also a signal peak                        
                    if ECG_bp_peaks[c] > THRESHOLDF1:                                            
                        SPKF = 0.125 * ECG_bp[c] + 0.875 * SPKF 
                        ECG_bp_signal_peaks.append(ECG_bp_peaks[c])                             
                    else:
                        ECG_bp_noise_peaks.append(ECG_bp_peaks[c])
                        NPKF = 0.125 * ECG_bp[c] + 0.875 * NPKF                   
                                        
            elif ECG_movavg[peak] > THRESHOLDI1 and ECG_movavg[peak] < THRESHOLDI2:
                #update noise thresholds
                NPKI = 0.125 * ECG_movavg[peak]  + 0.875 * NPKI  
                NPKF = 0.125 * ECG_bp[c] + 0.875 * NPKF
                    
            elif ECG_movavg[peak] < THRESHOLDI1:
                #update noise thresholds
                noise_peaks.append(peak)
                NPKI = 0.125 * ECG_movavg[peak]  + 0.875 * NPKI            
                ECG_bp_noise_peaks.append(ECG_bp_peaks[c])                       
                NPKF = 0.125 * ECG_bp[c] + 0.875 * NPKF
        else:
            if ECG_movavg[peak] >= THRESHOLDI1: #first peak is a signal peak
                IWF_signal_peaks.append(peak) 
                #update signal  thresholds
                SPKI = 0.125 * ECG_movavg[peak]  + 0.875 * SPKI
                if ECG_bp_peaks[c] > THRESHOLDF1:                                            
                    SPKF = 0.125 * ECG_bp[c] + 0.875 * SPKF 
                    ECG_bp_signal_peaks.append(ECG_bp_peaks[c])                             
                else:
                    ECG_bp_noise_peaks.append(ECG_bp_peaks[c])
                    NPKF = 0.125 * ECG_bp[c] + 0.875 * NPKF                                    
                
            elif ECG_movavg[peak] > THRESHOLDI2 and ECG_movavg[peak] < THRESHOLDI1:
                #update noise thresholds
                NPKI = 0.125 * ECG_movavg[peak]  + 0.875 * NPKI  
                NPKF = 0.125 * ECG[c] + 0.875 * NPKF
                                    
            elif ECG_movavg[peak] < THRESHOLDI2:
                #update noise thresholds
                noise_peaks.append(peak)
                NPKI = 0.125 * ECG_movavg[peak]  + 0.875 * NPKI            
                ECG_bp_noise_peaks.append(ECG_bp_peaks[c])                       
                NPKF = 0.125 * ECG_bp[c] + 0.875 * NPKF       
            
                    
        #reset 
        T_wave_found = 0                
            
        #update thresholds
        THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI)
        THRESHOLDI2 = 0.5 * THRESHOLDI1 
            
        THRESHOLDF1 = NPKF + 0.25 * (SPKF - NPKF)
        THRESHOLDF2 = 0.5 * THRESHOLDF1
    
    #adjust for filter delays
    ECG_R_locs = [int(i - delay) for i in ECG_bp_signal_peaks]
    ECG_R_locs = np.unique(ECG_R_locs)
    
    #neighborhood search in raw ECG signal for increase accuracy of R peak detection    
    for i in ECG_R_locs:
        ECG = np.array(ECG)
        searchInterval = int(np.round(0.02 * fs))
        searchIndices = np.arange(i - searchInterval, i + searchInterval + 1, 1)
        searchIndices = [i.item() for i in searchIndices] #convert to native Python int
        final_R_locs.append(np.where(ECG[searchIndices] == max(ECG[searchIndices]))[0][0] + searchIndices[0])
    
    #plot ECG signal with R peaks marked
    if plot == 1:
        samples = np.arange(0, len(ECG))
        plt.plot(samples, ECG, c = 'b')        
        plt.scatter(final_R_locs, ECG[final_R_locs], c = 'r', s = 30)
        plt.xlabel('Sample')
        plt.ylabel('ECG')
    else:
        pass
        
    return final_R_locs
                    
def findPeaks(ECG_movavg):
    """finds peaks in Integration Waveform by smoothing, locating zero crossings, and moving average amplitude thresholding"""
    #smoothing
    N = 15
    ECG_movavg_smooth = np.convolve(ECG_movavg, np.ones((N,)) / N, mode = 'same')    
    #signal derivative    
    sigDeriv = np.diff(ECG_movavg_smooth)     
    #find location of zero-crossings
    zeroCross = []
    for i,c in enumerate(np.arange(len(sigDeriv)-1)):
        if sigDeriv[i] > 0 and sigDeriv[i + 1] < 0:
            zeroCross.append(c)           
    
    return np.array(zeroCross) 


