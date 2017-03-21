# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 19:19:25 2017

@author: Sarah Pickus




"""

import numpy as np
from matplotlib import style
import matplotlib.pyplot as plt
style.use('ggplot')

def plotPoincare(RRints):
    
    #Implement Poincare plot
    ax1 = RRints[:-1]
    ax2 = RRints[1:]   
    plt.scatter(ax1, ax2, c = 'r', s = 15)
    plt.xlabel('RR_n (s)')
    plt.ylabel('RR_n+1 (s)')
    plt.show()

def eclipseFittingMethod(RRints):
    
    SDSD =  np.std(np.diff(RRints))
    SDRR = np.std(RRints)
    SD1 = (1 / np.sqrt(2)) * SDSD #measures the width of poincare cloud, indicating level of short term HRV
    SD2 = np.sqrt((2 * SDRR ** 2) - (0.5 * SDSD ** 2)) #measures the length of the poincare cloud, denotes long term HRV
 
    return {'SD1': SD1, 'SD2': SD2}
    
def asymmetryMethod(RRints):
    
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
    
    ax1 = RRints[:-1]
    ax2 = RRints[1:]
    mu_rr = np.mean(RRints)
    r_rr = np.mean((ax1 - mu_rr) * (ax2 - mu_rr)) / (np.sqrt(np.mean((ax1 - mu_rr) ** 2) * np.mean((ax2 - mu_rr) ** 2)))

    return r_rr
    
def plotRRintHist(RRints):
    
    plt.hist(RRints, bins = 'auto')
    plt.xlabel('RR Interval')
    plt.ylabel('Number of RR Intervals')
    plt.title('RR Interval Histogram')
    plt.show()

def plotWidthHist(RRints):
    
    ax1 = RRints[:-1]
    ax2 = RRints[1:]
    x1 = (np.cos(np.pi / 4) * ax1) - (np.sin(np.pi / 4) * ax2)
    plt.hist(x1, bins = 'auto')
    plt.title('Width (Delta-RR Interval) Histogram')
    plt.show()
    
def plotLengthHist(RRints):
    
    ax1 = RRints[:-1]
    ax2 = RRints[1:]
    x2 = (np.sin(np.pi / 4) * ax1) + (np.cos(np.pi / 4) * ax2)
    plt.hist(x2, bins = 'auto')
    plt.title('Length Histogram')
    plt.show()












