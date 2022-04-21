# Off-line and on-line Bayesian filtering

In sequential settings, Bayesian parameter estimation can be performed either in an off-line (batch) or an on-line (recursive) framework. This repository contains Matlab implementations of different Bayesian filters, whose goal is the quantification of the full posterior uncertainty of time-invariant model parameters in long-term monitoring settings. These are the following:
- Particle filter with Gaussian mixture resampling (PFGM)
- Particle filter with Gaussian mixture resampling and likelihood tempering (tPFGM)
- Independent Metropolis Hastings with Gaussian mixture proposal (IMH-GM)
- IMH-GM-based iterated batch importance sampling (IBIS) filter
- IMH-GM-based iterated batch importance sampling with likelihood tempering (tIBIS) filter
- IMH-GM-based sequential Monte Carlo filter

## Documentation

The following paper serves as documentation 

## Citation

Antonios Kamariotis, Luca Sardi, Iason Papaioannou, Eleni Chatzi, Daniel Straub (2022), On off-line and on-line Bayesian filtering for uncertainty quantification of structural deterioration, arXiv preprint.

## Notes

- Use of ERADist, ERANataf
