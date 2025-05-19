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
```
_______________________________________________________________________________________________________________________________________________________________________________

This code is employed for the estimation of time-varying (TV) linear and nonlinear systems modeled via Laguerre–Volterra expansions, with the aim of accurately recovering their evolving kernel structure using adaptive identification algorithms.
A synthetic system with mixed-mode variations—including both smooth (sinusoidal) and abrupt changes in the kernel coefficients—is simulated using CreateSimulation.m (a script that generates synthetic data and stores ground truth kernels (Realkernels)), with Gaussian white noise input and additive measurement noise in the output. The system itself is modeled using Laguerre basis expansions to represent first- or second-order Volterra kernels, whose time-varying coefficients capture the system's TV dynamics. Recursive algorithms are used to estimate these coefficients: RLS with constant forgetting (RLSC), RLS with adaptive forgetting (RLSA), the Kalman Filter (KF), and the Adaptive Kalman Filter (KFA). Each estimator operates on the regressor matrix constructed by GetV.m, with optional backward smoothing via kalsmooth.m or rlsmooth.m to reduce estimate variance. Optimal hyperparameters—such as model memory, forgetting factors, and noise covariances—are tuned using a Genetic Algorithm (GA) that minimizes either the Bayesian (BIC) or Akaike (AIC) information criterion. The RUNTHIS.m script controls the entire pipeline: generating the synthetic dataset, optimizing model parameters, applying the selected estimator, and visualizing both the predicted outputs and the true vs. estimated kernel surfaces. For more details see below.

________________________________________________________________________________________________________________________________________________________________________________
# References (please cite the following work, thank you!)

Kostoglou, Kyriaki, Ronald Schondorf, and Georgios D. Mitsis. "Modeling of multiple-input, time-varying systems with recursively estimated basis expansions." Signal Processing 155 (2019): 287-300.


________________________________________________________________________________________________________________________________________________________________________________


📦 Simulation and Kernel Generation

    CreateSimulation.m
    Simulates a time-varying nonlinear system with first- and optionally second-order Volterra kernels. Kernels vary over time using sinusoidal modulations with different     
    amplitudes/frequencies for each half of the signal. Adds noise to the output and saves the input, output, and true kernel.

    CREATEKERN.m
    Generates 1st- and optionally 2nd-order Volterra kernels using provided Laguerre basis matrices and coefficient vectors.

    Getlag.m
    Constructs Laguerre basis functions (lags) for given system memory length L, decay parameter a, and sample period T. Used to model kernel memory efficiently.

    GetV.m
    Builds the time-varying regressor matrix V from the input signal using Laguerre expansions. Includes terms for 1st- and optionally 2nd-order kernel modeling (Volterra 
    products).


🧠 Estimation Algorithms

Each of the following functions estimates the time-varying kernel coefficients from the input-output data using different adaptive methods. They all return a scalar cost J used in GA optimization.

    RLSC.m — Recursive Least Squares with constant forgetting factor

    RLSA.m — RLS with Adaptive forgetting (for tracking nonstationary systems)

    KF.m — Kalman Filter with fixed process and measurement noise

    KFA.m — Adaptive Kalman Filter where process noise variance updates over time

Each estimator is configured with parameters like L, a, lambda, R1, R2, etc., depending on the method.


🧬 Optimization

    RUNTHIS.m
    Main script: generates simulation, selects estimator method, runs Genetic Algorithm (GA) to optimize estimator hyperparameters, evaluates performance, and produces plots.     
    Controls the entire workflow.

    makesymcoef.m
    Symmetrizes the second-order kernel coefficients estimated by the recursive methods, ensuring proper symmetric structure (needed for Volterra kernels of order 2).


🧹 Smoothing

    kalsmooth.m
    Applies Kalman-style backward smoothing (Rauch–Tung–Striebel) to KF-based coefficient trajectories after forward estimation.

    rlsmooth.m
    Applies exponential backward smoothing to RLS-based coefficient trajectories using time-varying forgetting factors (useful to reduce high-frequency fluctuations in estimates).
