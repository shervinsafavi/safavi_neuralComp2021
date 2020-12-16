% this script plot figure 1 of the paper

%% 
figure
set(gcf,'defaultLineLineWidth',2.8)

% set the appearance of the figures
vc.f1.w = 16.3072;
vc.f1.h = 8.2176;
fh = gcf;
fh.Color        = [1 1 1]
fh.Units        = 'centimeters';
fh.Position(3)  = vc.f1.w;
fh.Position(4)  = vc.f1.h;

%% assign parameters 
% fix the seed of the random number generator to get consistent figure
rng(2);

% assign visualization parameters 
run(fullfile('vizConventions.m'))

% assign simulation parameters 
T      = 5;          % simulation length
dt     = .001;       % simulation step
t      = 0 : dt : T; % time samples
rate   = 8;          % event rate
nTrial = 1;          % number of trials in each simulation
nSim   = 1;          % number of simulation 
kappa  = .5;         % coupling strength

%% compute PLVs
clear PLV
clear nSpk

for ksim = 1 : nSim
    PLV(ksim) = 0;
    nSpk(ksim) = 0;
    for kTrial = 1:nTrial
        spkTimes = rand(length(t), 1) < (rate*exp(kappa * cos(2*pi*t')) * dt);
        % allspt : all spike train 
        allspt(kTrial, :) = spkTimes;
        PLV(ksim) = PLV(ksim) + sum(exp(1i*2*pi*t) .* spkTimes');
        nSpk(ksim) = nSpk(ksim)+sum(spkTimes);
    end
end

PLV = PLV ./ nSpk;


%% plot
plotRange = (1 : 3000) + 1000;

oscLambda = rate*exp(kappa*cos(2*pi*t));

hold all

s1 = plot(plotRange, 2*spkTimes(plotRange) - 7)
s2 = plot(plotRange, oscLambda(plotRange))
s3 = plot(plotRange, cumsum(spkTimes(plotRange)))
s4 = plot(plotRange, cumsum(spkTimes(plotRange))' - dt*cumsum(oscLambda(plotRange)), 'color',vc.f0.mc)

yline(0, 'linestyle','--', 'linewidth', 1, 'color',vc.f0.mc)

box off
axis off
ylim([-10 30])
legend([s1, s2, s3, s4], 'dN(t)','\lambda(t)','N(t)','M(t)', 'location', 'northwest')
set(gca, 'fontsize', 14)
xlabel('Time')


