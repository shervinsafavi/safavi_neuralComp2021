function [lfpLikeSig, spikeTrains, varargout] = smlt_vonMisesCouplingModel_wMultPop(globalDynamicsParams, spikeTrainParams, couplingParams, signalParams)
% similar to smlt_vonMisesCouplingMocel, but also LFP is consist of multiple oscillatory component instead of one.
%
% ------
% Input:
% 1     globalDynamicsParams: structure
%       characterize the global dynamics and should include the following fields:
% 1a        .oscFreq: scalar
%            oscilation frequency (in hertz) of locked spikes and LFPs
% 1b        .lfpPhaseNoise_kappa: scalar
%            TEXT
% 1c        .lockingStrength_kappa: scalar, vector or a matrix
%            parameter determining the concentration 
%            (similar to what we have in a von Mises distribution)
% (1d)      .lockingPhase: scalar, vector or a matrix
%            phase (in radian) where spikes are locked in    
% 1     spikeTrainParams: structure
%       contain features of the spike train(s) and should include the following fields:
% 2a        .avefiringRate: scalar, vector or a matrix
%            average spiking rate over time in the transient epochs                  
% 3     signalParams: structure
%       contain signal params and should include the following fields:
% 3a        .signalLength: scalar
%            length of your LFP time series and spike train in second
% 3b        .SF: scalar
%            sampling frequency
% 3c        .nTr: scalar
%            number of trials
% 3d        .nCh: scalar
%            number of channels
% 3e        .nUnit: scalar
%            number of spiking units 
% 
% Output:
% 1     lfpLikeSig: 
% 2     spikeTrains: binnary matrix 
%
% ------
% potential improvments:
% (1) TEXT
% ------
% Code Info:
%   creation: 2018-07-10 by SS (codes@shervinsafavi.org)
%   modification:
%       $ YYYY-MM-DD TEXT 
% ------
% see also

    %% some basic parameters
    signalParams.dt         = 1/signalParams.SF;               % length of time bin (in s)
    signalParams.t          = 0 : signalParams.dt : signalParams.signalLength - signalParams.dt;
    signalParams.nSample    = numel(signalParams.t);

    %% generate sustained (un-)coupled-spike

    for ifreq = 1 : numel(globalDynamicsParams.oscFreq)
        % PL = phase-locked
        PLspikeTrainParams.avefiringRate = spikeTrainParams.avefiringRate(:, ifreq); 
        PLspikeTrainParams.kappa        = couplingParams.lockingStrength_kappa(:,ifreq);
        PLspikeTrainParams.lockingFreq  = globalDynamicsParams.oscFreq(ifreq);
        PLspikeTrainParams.lockingPhase = couplingParams.lockingPhase;    
        [PLspikeTrains_allFreqSep(:,:,:, ifreq)]  = ...
            gnrt_phaseLockedSpikeTrains(PLspikeTrainParams, signalParams);
    end

    PLspikeTrains  = logical(sum(PLspikeTrains_allFreqSep, 4));
    spikeTrains_binary = PLspikeTrains;

    % convert the spike trains sparse matrices
    spikeTrains{signalParams.nTr} = [];
    for iTr = 1 : signalParams.nTr
        spikeTrains{iTr} = sparse(spikeTrains_binary(:,:, iTr));
    end

    %% generate continues osc LFP
    % probably it wont be needed (except the size of which is used later)
    noiseLFP = ...
        globalDynamicsParams.whiteNoise_sigma ...
        * randn(signalParams.nCh, signalParams.nSample, signalParams.nTr);

    nOscComp = numel(globalDynamicsParams.oscComps);
    amps = globalDynamicsParams.oscComps;
    oscWeights = amps / sum(amps);
    oscFreqs = globalDynamicsParams.oscFreq;


    for iOscComp = 1 : nOscComp
        cmplx_oscComps(iOscComp, :) = oscWeights(iOscComp) ...
            * exp(1i * (2*pi* oscFreqs(iOscComp) * signalParams.t));
    end

    cmplx_osc_raw = globalDynamicsParams.mixingMatrix * cmplx_oscComps;

    % oscWeights are useless, given we have mixing matrix here

    cmplx_noisyOsc = repmat(cmplx_osc_raw, 1, 1, signalParams.nTr) ...
        .* exp(1i * randvm(globalDynamicsParams.lfpPhaseNoise_kappa, size(noiseLFP))) ...
        + noiseLFP;

    lfpLikeSig = real(cmplx_noisyOsc);
end
