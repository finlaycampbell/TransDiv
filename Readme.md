# gensig: simulating genetic signatures of transmission events


### Introduction
Pathogen genetic sequence data can inform us of likely transmission trees if
genetic variation is introduced at a rate high enough to resolve individual
tranmsission pairs. We quantify the informativeness of genetic sequence data in
informing likely transmission pairs by calculating the 'genetic signature',
defined as the average number of mutations seperating a given transmission
pair.

Here we describe the genetic signatures of various viral and bacterial pathogens
by simulating sequence evolution in the context of a simple transmission model,
drawing values for R0, mutation rates, genome lengths, generation times and the
degree of heterogeneity in infectiousness from the literature. By counting the
number of mutations seperating individual transmission pairs, we explicitly
describe the distribution of values for the genetic signature in realistic
outbreak scenarios. This will inform us of epidemic situations for which genetic
sequence data should be considered in inferring transmission trees and informing
the implementation of infection control measures.

<br>

### Results
We simulated outbreaks and pathogen sequence evolution using the following
parameters:


|Pathogen      |  R0| Mutation rate<br>(base<sup>-1</sup> day<sup>-1</sup>)| Genome length| Mean generation time (days)| SD generation time (days)| Expected genetic signature|
|:-------------|---:|-----------------------------------------------------:|-------------:|---------------------------:|-------------------------:|--------------------------:|
|SARS          | 3.5|                                               7.6e-06|         29750|                         8.4|                       3.8|                       1.90|
|Klebsiella    | 2.0|                                               0.0e+00|       5305677|                        62.7|                      24.0|                       1.42|
|Ebola         | 5.0|                                               2.3e-06|         18058|                        15.3|                       9.3|                       0.63|
|MERS          | 3.5|                                               1.7e-06|         30115|                        10.7|                       6.0|                       0.54|
|Influenza     | 1.3|                                               7.9e-06|         13155|                         3.0|                       1.5|                       0.31|
|MRSA          | 1.3|                                               0.0e+00|       2842618|                        15.6|                      10.0|                       0.11|
|Streptococcus | 2.0|                                               0.0e+00|       2126652|                         6.6|                       1.8|                       0.05|

<br>
The distributions of realised generation times are given below.
![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)


### Contact
Finlay Campbell (f.campbell15@imperial.ac.uk) <br>
PhD Student <br>
Department of Infectious Disease Epidemiology <br>
Imperial College London <br>
