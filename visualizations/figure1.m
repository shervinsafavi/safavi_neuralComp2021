function figure1()
% figure1()
% this function plot the figure 1 of following paper:
% [From univariate to multivariate coupling between continuous signals and point processes: a mathematical framework, S.Safavi, N. K. Logothetis and M. Besserve. ArXiv 2020](https://arxiv.org/abs/2005.04034)


clf

%% assign parameters 

% fix the seed of the random number generator to get consistent figure
rng(2);


% assign visualization parameters (colors, etc)
vc = get_vizConventions();
set(gcf,'defaultLineLineWidth',2.8)
adjustFigAppearnce(vc.f1);

% assign parameters von Mises coupling model
vmCMparams.T      = 5;          % simulation length
vmCMparams.dt     = .001;       % simulation step
vmCMparams.t      = ...         % simulation time steps%
    0 : vmCMparams.dt : vmCMparams.T;      
vmCMparams.rate   = 8;          % event rate
vmCMparams.nTrial = 1;          % number of trials in each simulation
vmCMparams.nSim   = 1;          % number of simulation 
vmCMparams.kappa  = .5;         % coupling strength

%% simulate von Mises coupling model
[PLV, spkTimes] = smlt_vonMisesCouplingMocel(vmCMparams);

%%
% indices of time points used for plotting
tmpidx = (1000 : 4000);

% rate lambda
lambda = vmCMparams.rate * exp(vmCMparams.kappa * cos(2 * pi * vmCMparams.t));

hold all
s1 = plot(tmpidx, 2 * spkTimes(tmpidx) - 7)
% 2 and 7 here is just added to move it for visualization purposes
s2 = plot(tmpidx, lambda(tmpidx))
s3 = plot(tmpidx, cumsum(spkTimes(tmpidx)))
s4 = plot(tmpidx, cumsum(spkTimes(tmpidx)) - vmCMparams.dt*cumsum(lambda(tmpidx)), 'color',vc.f1.mc)
yline(0, 'linestyle','--', 'linewidth', 1, 'color', vc.f1.mc)

box off; axis off
ylim([-10 30])
legend([s1, s2, s3, s4], 'dN(t)','\lambda(t)','N(t)','M(t)', 'location', 'northwest')
set(gca, 'fontsize', 14)
xlabel('Time')

