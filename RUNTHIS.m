close all;clear all;

%Simulation parameters
N=1000;      %Number of data samples
fs=1;        %Sampling Rate
T=1/fs;
amp1=0.05;   %Amplitude of sinusoidal variations of the kernels (for the first half of the samples)
amp2=0.20;   %Amplitude of sinusoidal variations of the kernels (for the second half of the samples)
w1=2;        %Frequency of sinusoidal variations of the kernels (for the first half of the samples)
w2=10;       %Frequency of sinusoidal variations of the kernels (for the second half of the samples)
Q=1;         %Order of Nonlinearity (Q=1:1st-order,Q=2:2nd-order) 
Kerlen=100;  %Length of Volterra Kernel for representation (i.e. lags)
SNR=20;      %20dB SNR


%Optimization parameters for the GA
metric=1;   %1 for BIC / 2 for AIC
method=4;   %1 for RLSC / 2 for RLSA / 3 for KF / 4 for KFA
ignore=100; %ignore first 100 points of initialization

%Smoothing
smooth=1;   %1 to apply smoothing/0 for no smoothing


%------Create the simulation-----------------------------------------------
rng('shuffle')
inp1=randn(N,1);       %Input (Always a column vector)
CreateSimulation(inp1,Q,amp1,amp2,w1,w2,N,SNR,Kerlen)
ff=sprintf('SIM_w%d_w%d_amp%1.2f_amp%1.2f.mat',w1,w2,amp1,amp2);
load(ff)

%%%%%%%%%%%%%%%Identify TV System%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ga_opts = gaoptimset('TolFun',1e-12,'StallGenLimit',25,'Generations',100,'Display','iter');
%If you want to use parallel processing uncomment the following instead
% parpool(12)
% ga_opts = gaoptimset('TolFun',1e-12,'StallGenLimit',25,'Generations',50,'Display','iter','UseParallel','always');

if(method==1)
    nvars=4;
    LB=[1 0.1 0.1 0];
    UB=[10 0.9 1 Inf];
    h = @(X) RLSC(out,inp1,Q,X,metric,ignore,T);   %Use Matlab Coder and Create a mex function, which is faster. Once the mex function is created Use RLSC_mex instead of RLSC. The same applied for the other functions (i.e. RLSA, KF etc)
    rng('shuffle')
    [lam, err_ga] = ga(h, nvars,[],[],[],[],LB,UB,[],1,ga_opts);
elseif(method==2)
    nvars=6;
    LB=[1 0.1 0.1*ones(1,3) 0];
    UB=[10 0.9 ones(1,3) Inf];
    h = @(X) RLSA(out,inp1,Q,X,metric,ignore,T);
    rng('shuffle')
    [lam, err_ga] = ga(h, nvars,[],[],[],[],LB,UB,[],1,ga_opts);        
elseif(method==3)
    nvars=5;
    LB=[1 0.1 0 0 0];
    UB=[10 0.9 Inf Inf Inf];
    h = @(X)KF(out,inp1,Q,X,metric,ignore,T);
    rng('shuffle')
    [lam, err_ga] = ga(h, nvars,[],[],[],[],LB,UB,[],1,ga_opts);        
elseif(method==4)
    nvars=5;
    LB=[1 0.1 0.1 0 0];
    UB=[10 0.9 1 Inf Inf];
    h = @(X)KFA(out,inp1,Q,X,metric,ignore,T);
    rng('shuffle')
    [lam, err_ga] = ga(h, nvars,[],[],[],[],LB,UB,[],1,ga_opts);        
end

%----Print optimal GA solution---------------------------------------------
metr={'BIC','AIC'};
meth={'RLSC','RLSA','KF','KFA'};
if(method==1)
    fprintf('GA Optimal Solution based on %s and %s: L=%d a=%1.3f lambda=%1.3f P0=%1.3f', metr{metric},meth{method},lam)
elseif(method==2)
    fprintf('GA Optimal Solution based on %s and %s: L=%d a=%1.3f lambda_m_i_n=%1.3f lambda_m_a_x=%1.3f lambda_e=%1.3f P0=%1.3f', metr{metric},meth{method},lam)    
elseif(method==3)
    fprintf('GA Optimal Solution based on %s and %s: L=%d a=%1.3f R1=%1.3f R2=%1.3f P0=%1.3f', metr{metric},meth{method},lam) 
else
    fprintf('GA Optimal Solution based on %s and %s: L=%d a=%1.3f lambda_w=%1.3f R2=%1.3f P0=%1.3f', metr{metric},meth{method},lam)  
end

%----Evaluate GA solution--------------------------------------------------
results=TEST(out,inp1,Q,lam,ignore,T,method,Kerlen,smooth);
%-----Plots----------------------------------------------------------------
figure;plot(0:1/fs:(N-ignore)/fs,out(ignore:end));hold on;
plot(0:1/fs:(N-ignore)/fs,out_noisefree(ignore:end),'r')
plot(0:1/fs:(N-ignore)/fs,results.yhat,'g')
xlabel('n(sec)')
legend('Noisy output','Noise free output','Predicted output')

figure;colormap('jet');mesh(0:Kerlen-1,0:1/fs:(N-ignore)/fs,Realkernels.k1_1(ignore:end,:))
zlim=get(gca,'Zlim');
set(gca,'Zlim',[zlim(1)-0.5 zlim(2)+0.5])
ylabel('n(sec)')
xlabel('lags')
xlim([0 60])
title('Simulated 1st-order TV Laguerre Volterra Kernel') 
figure;colormap('jet');mesh(0:Kerlen-1,0:1/fs:(N-ignore)/fs,results.K1_1)
set(gca,'Zlim',[zlim(1)-0.5 zlim(2)+0.5])
ylabel('n(sec)')
xlabel('lags')
xlim([0 60])
title('Estimated 1st-order TV Laguerre Volterra Kernel') 
