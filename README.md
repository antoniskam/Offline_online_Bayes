# Off-line and on-line Bayesian filtering

In sequential settings, Bayesian parameter estimation can be performed either in an off-line (batch) or an on-line (recursive) framework. This repository contains Matlab implementations of different Bayesian filters, whose goal is the quantification of the full posterior uncertainty of time-invariant model parameters in long-term monitoring settings. These are the following:
- Particle filter with Gaussian mixture resampling (PFGM)
- Particle filter with Gaussian mixture resampling and likelihood tempering (tPFGM)
- Independent Metropolis Hastings with Gaussian mixture proposal (IMH-GM)
- IMH-GM-based iterated batch importance sampling (IBIS) filter
- IMH-GM-based iterated batch importance sampling with likelihood tempering (tIBIS) filter
- IMH-GM-based sequential Monte Carlo (SMC) filter

In this repository, the filters are applied on the following two numerical examples:
- A non-linear, non-Gaussian, low-dimensional probabilistic fatigue crack growth model that is updated with sequential crack monitoring measurements. The corresponding Matlab codes are found in the folder "cgm" of this repository.
- A linear, Gaussian, high-dimensional random field model of the spatially and temporally varying corrosion deterioration across a beam, which is updated with sequential measurements from sensors. The corresponding Matlab codes are found in the folder "hdRF" of this repository.

## Computational remarks
- The  algebraic  operations  in  all  presented  algorithms  are  implemented  in  the  logarithmic  scale,  which employs evaluations of the logarithm of the likelihood function and, hence, ensures computational stability.
- The  EM  step  for  fitting  the  Gaussian mixture model  is  performed  after  initially  transforming  the  prior  joint probability density function of the parameters to an underlying vector of independent standard normal random variables. This transformation is beneficial, since in the standard normal space the parameters are decorrelated,and that enhances the performance of the EM algorithm. This transformation is performed with the use of the ERANataf Matlab class, which is part of the ERADist Matlab package developed by the Engineering Risk Analysis (ERA) group at the Technical University of Munich. The ERADist software is freely available from this link (https://www.cee.ed.tum.de/era/software/). 

## Citation

Antonios Kamariotis, Luca Sardi, Iason Papaioannou, Eleni Chatzi, Daniel Straub (2022), On off-line and on-line Bayesian filtering for uncertainty quantification of structural deterioration, arXiv preprint.


