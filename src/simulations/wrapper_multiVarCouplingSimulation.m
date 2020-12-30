function [svdOut, varargout] = wrapper_multiVarCouplingSimulation(vmCMwMPparams, caseName)
    
    
    nUnitNum = numel(vmCMwMPparams.unitNums);
    % loop on different choice of units
    for iun = 1 : nUnitNum

        noc = numel(vmCMwMPparams.globalDynamicsParams(iun).oscFreq);

        %
        iCase = iun; % case based on number of units used

        
        signalParams(iCase) = vmCMwMPparams.signalParams;
        globalDynamicsParams = vmCMwMPparams.globalDynamicsParams(iCase);
        % iun = 1;
        % signalParams(iCase).nCh = 100;%unitNums(3); % use 100 channels
        %                               % signalParams.nUnit = unitNums(iun);
        signalParams(iCase).nUnit = vmCMwMPparams.unitNums(iCase);

        %%
        % cer = []; % clusterer: define the population structure
        % ncc = 1;  % we only play with one case of cell clusterting
        % nclu = globalDynamicsParams.nFreqComp; % number of clusters (because we deal with 6 freq)
        % cer = eye(nclu);
        % cer((3:6), (3:6)) = 0;

        nclu = globalDynamicsParams.nFreqComp;

        % iun = 1;
        % signalParams.nCh = unitNums(iun);
        % signalParams.nUnit = unitNums(iun);

        % coupling strength is specified later for the population
        couplingParams(iCase) = struct ...
            (...
                'lockingStrength_kappa', ... % concentration parameter $\kappa$ determining the strength of spike-LFP coupling
                NaN, ...
                'lockingPhase', ...          % locking phase of neurons (all neurons are lock-in a single phase)   
                2*pi * rand(signalParams(iCase).nUnit, 1) ... % rad  
                );

        % for icc = 1 : ncc % loop is not necessary here

        % mixing matrix
        baseLineCoef = vmCMwMPparams.mixingBaseLineCoef; % 0 correspond to no mixing
        globalDynamicsParams.mixingMatrix = ...
            kron(eye(noc), (1 - baseLineCoef)*ones(signalParams(iCase).nCh/noc,1)) + baseLineCoef;

        baseLineCoef = 0;
        spikeTrainParams.avefiringRate = kron(eye(noc), vmCMwMPparams.aveFR * (1 - baseLineCoef)*ones(signalParams(iCase).nUnit/noc,1)) + baseLineCoef;

        % ic = 7;
        % couplingStrength = couplingStrengths(ic);
        baseLineCoef = 0; % 0 correspond to no coupling across other freq
        nOscComp = numel(vmCMwMPparams.allOscFreq);
        tv = zeros(1, nOscComp);
        tv(1) = 1; % we are going to have 2 involved populations
        tv(end) = 1; % we are going to have 2 involved populations
        couplingPerPop = diag(tv);
        couplingParams(iCase).lockingStrength_kappa = kron(couplingPerPop, vmCMwMPparams.couplingStrength * (1 - baseLineCoef)*ones(signalParams(iCase).nUnit/ noc, 1)) + baseLineCoef;



        % icc = 1; % usless index, should be removed

        % ~ allocate 
        tmpSvdOut{1, vmCMwMPparams.nRel} = [];

        % for now
        statTestInfo = [];
        sameElecCheckInfo = [];
        
        parfor iRel = 1 : vmCMwMPparams.nRel

            %% simulation 
            [lfpLikeSig, spikeTrains] = ...
                smlt_vonMisesCouplingModel_wMultPop(...
                    globalDynamicsParams, spikeTrainParams, couplingParams(iCase), ...
                    signalParams(iCase));

            %% filter LFP 
            [~, analLfp] = pp_filt_recenter(lfpLikeSig, ...
                                            vmCMwMPparams.freqCenter, vmCMwMPparams.halfFilterWidth, ...
                                            signalParams(iCase).SF, ...
                                            vmCMwMPparams.filterOrder, vmCMwMPparams.nIteCent);

            %% GPLA >> need to be replaced with simpler routine

            % [~,~,~,~, allStats, tmpSvdOut{1, iRel}] = ...
            %     tngpla(spikeTrains, analLfp, [], [], [], [] , [], ...
            %            statTestInfo, 'all', sameElecCheckInfo, ...
            %            'nSpk-square-root', 1, 0);
            
            tmpSvdOut{1, iRel} = multVarCouplingAnalysis(spikeTrains, analLfp);

            % tmpSummaryStat(iun, ic, iRel) = allStats.gPLV_stats.nullHypoReject;
        end                                     %
        svdOut(iCase).(caseName) = tmpSvdOut;
        % svdOut{iCase} = tmpSvdOut;
        
        % fn = fieldnames(tmpSvdOut);
        % for kf = 1 : numel(fn)
        %     svdOut(iCase).(caseName)(fn{kf}) = tmpSvdOut.(fn{kf})
        % end
        
    end
    
    vmCMwMPparams.signalParams = signalParams;
    varargout{1} = vmCMwMPparams;
    
end


function [recentLfpPh, varargout] = pp_filt_recenter(lfpLikeSig, freqCenter, halfFilterWidth, ...
                                                     SF, filterOrder, nIteCent)
    % this function do the filtering and phase recursive recentering of
    % LFP data

    signalParams.nTr = size(lfpLikeSig, 3);
    signalParams.nCh = size(lfpLikeSig, 1);

    %% Filter parameters
    % spike-lfp analysis
    % freqBand = globalDynamicsParams.oscFreq(iFreq) + [-7.5 +7.5];
    freqBand = freqCenter + [-halfFilterWidth +halfFilterWidth];
    wn = freqBand / (SF/2);
    % filterOrder = 2;
    [b, a] = butter(filterOrder, wn, 'bandpass');


    %% filter LFP
    for iTr = 1 : signalParams.nTr
        for iCh = 1 : signalParams.nCh
            tmpLfp = lfpLikeSig(iCh,:, iTr);
            tmpLfp_filtered = filtfilt(b, a, tmpLfp);
            tmpAnlcLfp = hilbert(tmpLfp_filtered);
            
            analLfp(iCh,:, iTr) = tmpAnlcLfp;
            
            for kite = 1 :  nIteCent
                tmpAnlcLfp = tmpAnlcLfp./ abs(tmpAnlcLfp);
                tmpAnlcLfp = tmpAnlcLfp - mean(tmpAnlcLfp);
            end 
            
            recentLfpPh(iCh,:, iTr) = angle(tmpAnlcLfp);
            
        end
    end


    varargout{1} = analLfp;
end