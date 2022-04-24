# Off-line and on-line Bayesian filtering
In sequential settings, Bayesian parameter estimation can be performed either in an off-line (batch) or an on-line (recursive) framework. This repository contains Matlab implementations of different Bayesian filters, whose goal is the quantification of the full posterior uncertainty of time-invariant model parameters in long-term monitoring settings. These are the following:
- On-line Particle Filter with Gaussian Mixture resampling (PFGM).
- On-line Particle Filter with Gaussian Mixture resampling and likelihood tempering (tPFGM).
- On-line Iterated Batch Importance Sampling (IBIS) filter.
- On-line Iterated Batch Importance Sampling filter with likelihood tempering (tIBIS).
- Off-line Sequential Monte Carlo (SMC) filter, which applies likelihood tempering to sequentially arrive to a single posterior density of interest.
- The IBIS, tIBIS and SMC filters employ an MCMC move step via Independent Metropolis Hastings with Gaussian mixture proposal (IMH-GM).

In this repository, the filters are applied on the following two numerical examples:
- cgm: a non-linear, non-Gaussian, low-dimensional probabilistic fatigue crack growth model that is updated with sequential crack monitoring measurements. 
- hdRF: a linear, Gaussian, high-dimensional random field model of the spatially and temporally varying corrosion deterioration across a beam, which is updated with sequential measurements from sensors.

## Computational remark
- The  EM  step  for  fitting  the  Gaussian mixture model  is  performed  after  initially  transforming  the  prior  joint probability density function of the parameters to an underlying vector of independent standard normal random variables. This transformation is beneficial, since in the standard normal space the parameters are decorrelated,and that enhances the performance of the EM algorithm. This transformation is performed with the use of the ERANataf Matlab class, which is part of the ERADist Matlab package developed by the Engineering Risk Analysis (ERA) group at the Technical University of Munich. The ERADist software is freely available from this link (https://www.cee.ed.tum.de/era/software/). 

## Citation
This repository contains the algorithm implementations and numerical examples presented in this publication:

- Antonios Kamariotis, Luca Sardi, Iason Papaioannou, Eleni Chatzi, Daniel Straub (2022), On off-line and on-line Bayesian filtering for uncertainty quantification of structural deterioration, arXiv preprint.


