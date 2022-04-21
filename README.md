# Off-line and on-line Bayesian filtering

In sequential settings, Bayesian parameter estimation can be performed either in an off-line (batch) or an on-line (recursive) framework. This repository contains Matlab implementations of different Bayesian filters, whose goal is the quantification of the full posterior uncertainty of time-invariant model parameters in long-term monitoring settings. These are the following:
- Particle filter with Gaussian mixture resampling (PFGM)
- Particle filter with Gaussian mixture resampling and likelihood tempering (tPFGM)
- Independent Metropolis Hastings with Gaussian mixture proposal (IMH-GM)
- IMH-GM-based iterated batch importance sampling (IBIS) filter
- IMH-GM-based iterated batch importance sampling with likelihood tempering (tIBIS) filter
- IMH-GM-based sequential Monte Carlo filter

In this repository, the filters are applied on the following two numerical examples:
- A non-linear, non-Gaussian, low-dimensional probabilistic fatigue crack growth model that is updated with sequential crack monitoring measurements. The corresponding Matlab codes are found in the folder "cgm" of this repository.
- A linear, Gaussian, high-dimensional random field model of the spatially and temporally varying corrosion deterioration across a beam, which is updated with sequential measurements from sensors. The corresponding Matlab codes are found in the folder "hdRF" of this repository.

## Notes

Use of ERADist, ERANataf


## Citation

Antonios Kamariotis, Luca Sardi, Iason Papaioannou, Eleni Chatzi, Daniel Straub (2022), On off-line and on-line Bayesian filtering for uncertainty quantification of structural deterioration, arXiv preprint.


