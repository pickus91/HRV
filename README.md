# HRV

This is a Python module/tutorial for performing heart rate variability (HRV) analysis on electrocardiogram (ECG) time series. 

## Dependencies
* [NumPy](http://www.numpy.org/)
* [Matplotlib](http://matplotlib.org/)
* [SciPy](https://www.scipy.org/)

## Licence
This code is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## R-Peak Detection
### Pan-Tompkins

Implements the popular QRS complex detection algorithm introduced in [Pan, *et al* (1985)](https://www.researchgate.net/profile/Keesam_Jeong/publication/3728672_A_simple_real-time_QRS_detection_algorithm/links/54e829e10cf2f7aa4d4f64a9.pdf). The algorithm uses filtering, adaptive thresholding, and criteria based on human cardiac physiology to detect QRS complexes in the face of noise and quickly changing and diverse ECG morphologies. This function implements the Pan-Tompkins algorithm as it was originally published, along with two modifications which include additional filtering and eliminating potential QRS detections that occur within the refractory period. Since this algorithm is often used to find R-peak locations (and not just general QRS detection) for applications such as Heart Rate Variability (HRV) analysis, this function also performs a neighborhood search around the final QRS detection locations to find exact R-peak locations.

#### Code Example
```
    R_peak_locs = panTompkins(ECG, fs, plot = 1)

```

<div align = "center">
<img style="float: left;" src="https://github.com/pickus91/HRV/blob/master/figures/Original%20Signal.png"  height="350" width="425">
<img style="float: right;" src="https://github.com/pickus91/HRV/blob/master/figures/Final%20R%20Peak%20detection.png"  height="350" width="425">
</div>



## HRV Features

### Time Domain 
Time domain features are perhaps the simplist method computationally for performning HRV analysis. Following R-peak detection (see [PanTompkins](https://github.com/pickus91/HRV/blob/master/panTompkins.py)), ectopic beats are removed to produce a normal-to-normal (NN) interval series. Various statistical measures can then be extracted from the time series that provide insights into heart rate linear dynamics. 

#### Code Example

```
timeDomainFeats = timDomain(RR_interval_series)
```

| Label         | Description                                                       |
|:-------------:| :---------------------------------------------------------------- |
| ANN           | Average NN interval                                               | 
| SDNN          | Standard deviation of NN intervals                                |   
| SDSD          | Standard deviation of successive NN intervals                     | 
| NN50          | Number of successive NN intervals differing by more than 50 ms    |
| pNN50         | Proportion of successive NN intervals differing by more than 50 ms|
| rMMSD         | Root mean square of successive NN intervals                       |
| MedianNN      | Median of NN intervals                                            |


### Frequency Domain 
Spectral analysis is a standard in HRV analysis. Features are extracted from the power spectral density (PSD) of the NN interval time series. The PSD within various frequency bands is estimated using [Welch's method](https://en.wikipedia.org/wiki/Welch%27s_method). The most common frequency bands under investigation are the very low frequency (VLF) band, the low frequency (LF) band, and the high frequency (HF) band. The spectral power distributions across these bands has been shown to provide insight into modulations in autnomic nervous system (ANS) activity.

#### Code Example

```
freqDomainFeats = frequencyDomain(RR_interval_series)
```

  | Label         | Description                                                      |
  |:-------------:|:---------------------------------------------------------------- |
  | VLF Power     | Log of normalized spectral power between 0.003 Hz and 0.04 Hz    | 
  | LF Power      | Log of normalized spectral power between 0.04 Hz and 0.15 Hz     |   
  | HF Power      | Log of normalized spectral power between 0.15 Hz and 0.4 Hz      | 
  | LF/HF Ratio   | Ratio between LF and HF spectral power                           |


<div align = "center">
<img src="https://github.com/pickus91/HRV/blob/master/figures/frequencyDomain.png" align = "center" height="350" width="450"> 
</div>


The ```frequencyDomain``` function also has the capability to compute the frequency domain features listed in the table above using the spectral boundary adaptation method described in [X. Long, *et al* (2014)](https://pure.tue.nl/ws/files/3855045/5718586174038081.pdf). The idea behind using adapted spectral bands is to account for the time varying behavior of these bands (e.g. whilst sleeping) and to potentially reduce the within – and between – subject frequency domain feature variabilities. This method works by adapting the bands according to the peak frequencies found within the traditional spectral bands and setting their lower and upper bounds to values based off specified bandwidths.

<div>
<ul>        
<img style="float: left;" src="https://github.com/pickus91/HRV/blob/master/figures/frequencyDomain_traditional_bands.png"  height="350" width="400">
<img style="float: right;" src="https://github.com/pickus91/HRV/blob/master/figures/frequencyDomain_adapted_bands.png"  height="350" width="400">
 </ul>
</div>

### Poincare 

Poincare plots are an important visualization technique for quantifying the non-linear characteristics of the RR interval time series. The geometrical descriptors that can be extracted from the poincare plot with this package have been shown to provide insights into short and long term HRV trends. This includes parameters derived from the ellipse fitting method described in [Brennan, *et al* (2001)](http://ieeexplore.ieee.org/abstract/document/1018984/) and the heart rate asymmetry (HRA) method that quantifies accelerations/decelerations in heart rate introduced in [Guzik, Przemyslaw, *et al* (2006)](https://www.researchgate.net/profile/Przemyslaw_Guzik/publication/6734042_Heart_rate_asymmetry_by_Poincare_plots_of_RR_intervals/links/00463516712a5287a9000000/Heart-rate-asymmetry-by-Poincare-plots-of-RR-intervals.pdf).

<div align = "center">
<img src="https://github.com/pickus91/HRV/blob/master/figures/PoincarePlot.png" align="center" height="350" width="450">
</div>

### Multifractal Detrended Fluctuation Analysis (MF-DFA)
Multi-fractal detrended fluctuation analysis (MF-DFA) introduced in [Kantelhardt, *et al*](https://arxiv.org/pdf/physics/0202070.pdf).
MF-DFA is based on the standard detrended fluctuation analysis (DFA) introduced by [Peng, *et al* (1995)](http://havlin.biu.ac.il/PS/Quantification%20of%20scaling%20exponents%20and%20crossover%20phenomena%20in%20nonstationary%20heartbeat%20time%20series.pdf). 

#### Code Example
```
from DFA import scalingExponent
scaleExp = scalingExponent(timeSeries, lowerScaleLimit = 10, upperScaleLimit = 40, scaleDensity = 30, m = 1, q = 2, plot = 1)

```


The algorithm involves dividing the integrated RR interval time series into non-overlapping segments of equal length, *s*. 
<div align = "center">
<img src="https://github.com/pickus91/HRV/blob/master/figures/Original%20RR%20Series.png" height="300" width="350">
<img src="https://github.com/pickus91/HRV/blob/master/figures/Step%201%20-%20Integrated%20Time%20Series.png" height="300" width="350">
</div>
A polynomial is then fitted to each non-overlapping segment. Linear, quadratic, cubic, and/or other higher polynomials may be used in the fitting procedure. These DFA orders differ in their ability to eliminate various types of trends from the time series. Thus, an estimation of the type of time series trend can be captured via comparison of the different DFA orders, as seen below. 

<div>
        <ul>
            <img src="https://github.com/pickus91/HRV/blob/master/figures/Step%203%20-%20Calculate%20local%20trend%20-%20m%20%3D%201.png" alt="" class="first" height = "250" width = "275">
            <img src="https://github.com/pickus91/HRV/blob/master/figures/Step%203%20-%20Calculate%20local%20trend%20-%20m%20%3D%202.png" alt="" class="" height = "250" width = "275">
            <img src="https://github.com/pickus91/HRV/blob/master/figures/Step%203%20-%20Calculate%20local%20trend%20-%20m%20%3D%203.png" alt="" class="" height = "250" width = "275">              
        </ul>
<div>

Following the polynomial fit, the average of all the segments are obtained via:
<div align = "center">
<img src = "https://github.com/pickus91/HRV/blob/master/figures/Fluctuation%20Coefficient.PNG" align="center" height = "75" width = "200">
</div>

where *q* is the order of the fluctuation coefficient. When *q* = 2, the MF-DFA procedure simplifies to standard DFA. It may be of the user’s interests to explore how the *q*-dependent fluctuation coefficients *F<sub>q</sub>(s)* depend on scale *s* for various values of *q*. The above procedure is repeated over various scales to provide a relationship between the qth order fluctuation coefficient and scale. The final step of the MF-DFA procedure is determining the scaling behavior of the various fluctuation coefficients by generating a log-log plot of *F<sub>q</sub>(s)* versus *s*. In general, *F<sub>q</sub>(s)* increases with increases in scale, with a linear relationship on the double log plot indicating the presence of scaling. This behavior can be characterized by a scaling exponent α, which is the slope of the line of best fit describing the relationship between *F<sub>q</sub>(s)* and scale.

<div align = "center">
<img src = "https://github.com/pickus91/HRV/blob/master/figures/DFA%20Output.png" align = "center" height = "400" width = "500">
</div>

### Multiscale Entropy (MSE)

Implements the multiscale entropy (MSE) algorithm introduced in [Costa, *et al* (2002)](http://dbiom.org/files/publications/Peng_MultiscaleEntropyAnalysisComplexPhysiologicTimeSeries.pdf). MSE can be used as a tool for measuring the complexity of the RR tachogram. The MSE algorithm works by computing sample entropy at multiple scales (for various "coarse-grained" series), making it advantageious for analyzing features related to structure on scales other than the shortest one. Details on computing sample entropy are detailed in [Pincas, *et al* (1991)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC51218/pdf/pnas01056-0271.pdf). 

#### Code Example
```
import multiScaleEntropy as MSE
entropyMeasure = MSE.multiScaleEntropy(RR_interval_series,
                                       mseScales = np.arange(1,11),
                                       r = 0.2 * np.std(RRints),
                                       m = 2)
```

As seen below, applying the MSE algorithm to (uncorrelated) white noise results in a monotonically decreasing entropy measure, while its application to pink noise (or 1/f noise) reveals structure across multiple scales.

<div align = "center">
<img src = "https://github.com/pickus91/HRV/blob/master/figures/Noise_EntropyMeasures.png" align = "center" height = "400" width = "500">
</div>


