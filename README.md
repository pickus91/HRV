
# HRV

This is a Python module for performing heart rate variability (HRV) analysis on electrocardiogram (ECG) time series. 

## Dependencies

* [NumPy] (http://www.numpy.org/)
* [Matplotlib] (http://matplotlib.org/)
* [SciPy] (https://www.scipy.org/)

## R-Peak Detection
### Pan-Tompkins

Implements the popular QRS complex detection algorithm introduced in:

[Pan, Jiapu, and Willis J. Tompkins. "A real-time QRS detection algorithm." *IEEE transactions on biomedical engineering* 3 (1985): 230-236.](https://www.researchgate.net/profile/Keesam_Jeong/publication/3728672_A_simple_real-time_QRS_detection_algorithm/links/54e829e10cf2f7aa4d4f64a9.pdf)

The algorithm uses filtering, adaptive thresholding, and criteria based on human cardiac physiology to detect QRS complexes in the face of noise and quickly changing and diverse ECG morphologies. This function implements the Pan-Tompkins algorithm as it was originally published, along with two modifications which include additional filtering and eliminating potential QRS detections that occur
within the refractory period.Since this algorithm is often used to find R-peak locations (and not just general QRS detection) for applications such as Heart Rate Variability (HRV) analysis, this function also performs a neighborhood search
around the final QRS detection locations to find exact R-peak locations.

<img src="https://github.com/pickus91/HRV/blob/master/figures/Original%20Signal.png" align="center" height="350" width="450">
<img src="https://github.com/pickus91/HRV/blob/master/figures/Final%20R%20Peak%20detection_Init%20Phase.PNG" align="center" height="350" width="450">

## HRV Features

### Time Domain 

| Label         | Description                                                       |
|:-------------:| :---------------------------------------------------------------- |
| ANN           | Average NN interval                                               | 
| SDNN          | Standard deviation of NN intervals                                |   
| SDSD          | Standard deviation of successive NN intervals                     | 
| NN50          | Number of successive NN intervals differing by more than 50 ms    |
| pNN50         | Proportion of successive NN intervals differing by more than 50 ms|
| rMMSD         | Root mean square of successive NN intervals                       |
| MedianNN      | Median of NN intervals                                            |

### Frequency Domain 

  | Label         | Description                                                      |
  |:-------------:|:---------------------------------------------------------------- |
  | VLF Power     | Log of normalized spectral power between 0.003 Hz and 0.04 Hz    | 
  | LF Power      | Log of normalized spectral power between 0.04 Hz and 0.15 Hz     |   
  | HF Power      | Log of normalized spectral power between 0.15 Hz and 0.4 Hz      | 
  | LF/HF Ratio   | Ratio between LF and HF spectral power                           |
  
<img src="https://github.com/pickus91/HRV/blob/master/figures/frequencyDomain.png" align="center" height="350" width="450"> 

### Poincare 

Poincare plots are an important visualization technique for quantifying the non-linear characteristics of the RR interval time series. The geometrical descriptors that can be extracted from the poincare plot with this package have been shown to provide insights into short and long term HRV trends. This includes parameters derived from the ellipse fitting method described in [Brennan, *et al* (2001)](http://ieeexplore.ieee.org/abstract/document/1018984/) and the heart rate asymmetry (HRA) method that quantifies accelerations/decelerations in heart rate introduced in [Guzik, Przemyslaw, *et al* (2006)](https://www.researchgate.net/profile/Przemyslaw_Guzik/publication/6734042_Heart_rate_asymmetry_by_Poincare_plots_of_RR_intervals/links/00463516712a5287a9000000/Heart-rate-asymmetry-by-Poincare-plots-of-RR-intervals.pdf)

<img src="https://github.com/pickus91/HRV/blob/master/figures/PoincarePlot.png" align="center" height="350" width="450">

### Multifractal Detrended Fluctuation Analysis (MF-DFA)
Multi-fractal detrended fluctuation analysis (MF-DFA) introduced in:

[Kantelhardt, Jan W., et al. "Multifractal detrended fluctuation analysis of 
nonstationary time series." *Physica A: Statistical Mechanics and its Applications*
316.1 (2002): 87-114] (https://arxiv.org/pdf/physics/0202070.pdf)

MF-DFA is based on the standard detrended fluctuation analysis (DFA) introduced by [Peng, *et al* (1995)](http://havlin.biu.ac.il/PS/Quantification%20of%20scaling%20exponents%20and%20crossover%20phenomena%20in%20nonstationary%20heartbeat%20time%20series.pdf). The algorithm involves dividing the integrated RR interval time series into non-overlapping segments of equal length, *s*. 

<img src="https://github.com/pickus91/HRV/blob/master/figures/Original%20RR%20Series.png" align="center" height="300" width="350">

<img src="https://github.com/pickus91/HRV/blob/master/figures/Step%201%20-%20Integrated%20Time%20Series.png" align="center" height="300" width="350">

A polynomial is then fitted to each non-overlapping segment. Linear, quadratic, cubic, and/or other higher polynomials may be used in the fitting procedure. These DFA orders differ in their ability to eliminate various types of trends from the time series. Thus, an estimation of the type of time series trend can be captured via comparison of the different DFA orders, as seen below. 



          



<div>
        <ul>
            <img src="https://github.com/pickus91/HRV/blob/master/figures/Step%203%20-%20Calculate%20local%20trend%20-%20m%20%3D%201.png" alt="" class="first" height = "250" width = "300" align = "center">
            <img src="https://github.com/pickus91/HRV/blob/master/figures/Step%203%20-%20Calculate%20local%20trend%20-%20m%20%3D%202.png" alt="" class="" height = "250" width = "300" align = "center">
            <img src="https://github.com/pickus91/HRV/blob/master/figures/Step%203%20-%20Calculate%20local%20trend%20-%20m%20%3D%203.png" alt="" class="" height = "250" width = "300" align = "center" >              
        </ul>
<div>



Following the polynomial fit, the average of all the segments are obtained via:

<img src = "https://github.com/pickus91/HRV/blob/master/figures/Fluctuation%20Coefficient.PNG" align = "center" height = "75" width = "200">

where *q* is the order of the fluctuation coefficient. When *q* = 2, the MF-DFA procedure simplifies to standard DFA. It may be of the user’s interests to explore how the q-dependent fluctuation coefficients *F<sub>q</sub>(s)* depend on scale s for various values of *q*.

The above procedure is repeated over various scales to provide a relationship between the qth order fluctuation coefficient and scale. The final step of the MF-DFA procedure is determining the scaling behavior of the various fluctuation coefficients by generating a log-log plot of *F<sub>q</sub>(s)* versus *s*. In general, *F<sub>q</sub>(s)* increases with increases in scale, with a linear relationship on the double log plot indicating the presence of scaling. This behavior can be characterized by a scaling exponent α, which is the slope of the line of best fit describing the relationship between *F<sub>q</sub>(s)* and scale.

<img src = "https://github.com/pickus91/HRV/blob/master/figures/DFA%20Output.png" align = "center" height = "300" width = "350">

### Multiscale Entropy (MSE)

Costa, Madalena, Ary L. Goldberger, and C-K. Peng. "Multiscale entropy 
analysis of complex physiologic time series." *Physical review letters* 89.6 
(2002): 068102.
