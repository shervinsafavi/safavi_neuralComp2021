function [PLV, spikeTrains] = smlt_vonMisesCouplingMocel(vmCMparams)
% [PLV, spikeTrains] = smlt_vonMisesCouplingMocel(vmCMparams)
% This function simulate von Mises Coupling Mocel.
% See https://arxiv.org/abs/2005.04034 for more details.

I_0 = @(k) integral(@(x) exp(k * cos(x)), -pi, pi) / (2*pi);

for ksim = 1 : vmCMparams.nSim
    PLV(ksim) = 0;
    nSpk(ksim) = 0;
    for kTrial = 1 : vmCMparams.nTrial
        spkTimes = rand(length(vmCMparams.t), 1) ...
            < (vmCMparams.rate * exp(vmCMparams.kappa * cos(2 * pi * vmCMparams.t')) * vmCMparams.dt/ I_0(vmCMparams.kappa));
        % with division by vmCMparams.kappa, vmCMparams.rate is the average firing rate

        % 
        % spkTimes = rand(length(t),1)<(rate*exp(kappa*cos(2*pi*t'))*dt / I_0(kappa));
        % spikeTrains : all spike train 
        spikeTrains(kTrial, :) = spkTimes;
        PLV(ksim) = PLV(ksim) + sum(exp(1i * 2 * pi * vmCMparams.t) .* spkTimes');
        nSpk(ksim) = nSpk(ksim)+sum(spkTimes);
    end
end

PLV = PLV ./ nSpk;
