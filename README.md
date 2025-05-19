# TV-LV
Time-varying Laguerre Volterra models 
```
├── code
│   └── [matlab scripts]
│        RUNTHIS.m      -----> main function
│        CreateSimulation.m
│        CREATEKERN.m  
│        makesymcoef.m 
│        Getlag.m 
│        GetV.m 
│        RLSA.m 
│        RLSC.m 
│        rlsmooth.m 
│        KF.m 
│        KFA.m 
│        kalsmooth.m 
│        TEST.m 
│
├── tutorial
     └── tutorial.pdf
```
_______________________________________________________________________________________________________________________________________________________________________________

This code is employed for the estimation of time-varying (TV) linear and nonlinear systems modeled via Laguerre–Volterra expansions, with the aim of accurately recovering their evolving kernel structure using adaptive identification algorithms.
A synthetic system with mixed-mode variations—including both smooth (sinusoidal) and abrupt changes in the kernel coefficients—is simulated using CreateSimulation.m (a script that generates synthetic data and stores ground truth kernels (Realkernels)), with Gaussian white noise input and additive measurement noise in the output. The system itself is modeled using Laguerre basis expansions to represent first- or second-order Volterra kernels, whose time-varying coefficients capture the system's TV dynamics. Recursive algorithms are used to estimate these coefficients: RLS with constant forgetting (RLSC), RLS with adaptive forgetting (RLSA), the Kalman Filter (KF), and the Adaptive Kalman Filter (KFA). Each estimator operates on the regressor matrix constructed by GetV.m, with optional backward smoothing via kalsmooth.m or rlsmooth.m to reduce estimate variance. Optimal hyperparameters—such as model memory, forgetting factors, and noise covariances—are tuned using a Genetic Algorithm (GA) that minimizes either the Bayesian (BIC) or Akaike (AIC) information criterion. The RUNTHIS.m script controls the entire pipeline: generating the synthetic dataset, optimizing model parameters, applying the selected estimator, and visualizing both the predicted outputs and the true vs. estimated kernel surfaces.

________________________________________________________________________________________________________________________________________________________________________________
# References (please cite the following work, thank you!)

Kostoglou, Kyriaki, Ronald Schondorf, and Georgios D. Mitsis. "Modeling of multiple-input, time-varying systems with recursively estimated basis expansions." Signal Processing 155 (2019): 287-300.
