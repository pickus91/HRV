# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 13:59:24 2017

@author: picku
"""

import numpy as np

def timeDomain(NN):
    
    L = len(NN)    
    ANN = np.mean(NN)
    SDNN = np.std(NN)
    SDSD = np.std(np.diff(NN))    
    NN50 = len(np.where(np.diff(NN) > 0.05)[0])    
    pNN50 = NN50/L    
    NN20 = len(np.where(np.diff(NN) > 0.02)[0])
    pNN20 = NN20/L
    rMSSD = np.sqrt((1/L) * sum(np.diff(NN) ** 2))        
    MedianNN = np.median(NN)
    
    timeDomainFeats = {'ANN': ANN, 'SDNN': SDNN,
                       'SDSD': SDSD, 'NN50': NN50,
                       'pNN50': pNN50, 'NN20': NN20,
                       'pNN20': pNN20, 'rMSSD': rMSSD,
                       'MedianNN':MedianNN}
                       
    return timeDomainFeats