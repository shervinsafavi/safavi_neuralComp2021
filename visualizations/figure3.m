function figure3()
% figure3()
% this function plot the figure 3 of following paper:
% [From univariate to multivariate coupling between continuous signals and point processes: a mathematical framework, S.Safavi, N. K. Logothetis and M. Besserve. ArXiv 2020](https://arxiv.org/abs/2005.04034)

    clf

    %% assign parameters

    % fix the seed of the random number generator to get consistent figure
    rng(2);

    % assign visualization parameters 
    vc = get_vizConventions();
    adjustFigAppearnce(vc.f3);

    nR = vc.f3.nR;
    nC = vc.f3.nC;

    mlw = 2;
    ms = 10;
    fscp = 10;

    % assign parameters for generating homogeneous Poisson spikes and accompanying oscillation 
    % hPSandOscParams: homogeneous Poisson spikes and [accompanying] oscillation parameters 
    
    hPSandOscParams.T      = NaN;        % simulation length - is defined later 
    hPSandOscParams.dt     = .001;       % simulation step
    hPSandOscParams.rate   = 30;         % event rate
    hPSandOscParams.nTrial = 10;         % number of trials in each simulation
    hPSandOscParams.nSim   = 500;        % number of simulation 
    hPSandOscParams.f   = 1;             % frequency of accompanying oscillation 
    
    gamma_T_vals = 0.25 : 0.25 : 5;
    % gamma_T_vals reflect different choice for the length of the signal
    % length of the spike train are gamma_T * T 
    % see corollary 3 and Eq 18 of https://arxiv.org/abs/2005.04034

    %% simulation of homogeneous Poisson spike train and compute PLVs with perfectly linear phase

    % loop on different choice of signal length
    for kAlphaValIndex = 1 : numel(gamma_T_vals)
        gamma_T = gamma_T_vals(kAlphaValIndex);
        hPSandOscParams.gamma_T = gamma_T;

        % Inside "cmpt_PLV_wPoissonSpkAndLinPhase", gamma_T is used to calculate T
        for k = 1 : hPSandOscParams.nSim
            PLV(k, kAlphaValIndex) = cmpt_PLV_wPoissonSpkAndLinPhase(hPSandOscParams);
            PLV_groundTruth(kAlphaValIndex) = (2*pi*gamma_T*1i)^-1 * (exp(2*pi*gamma_T*1i) - 1);
            % ground trugh is calculated based on Eq 18 of https://arxiv.org/abs/2005.04034
        end
    end

    %% plotting PLVs
    subplot2d(nR,nC, 2, [1 2 3])
    boxplot(abs(PLV), gamma_T_vals', 'OutlierSize',1)
    hold on
    theoFH = plot(abs(PLV_groundTruth), ':', 'linewidth', 1.5, 'color',vc.gtc)

    k = 1;
    intPerLoc = 4 : 4 : numel(gamma_T_vals);
    intPerFH = xline(intPerLoc(k), 'b--');
    for k = 2 : numel(intPerLoc)%max(gamma_T_vals)
        xline(intPerLoc(k), 'b--');
    end

    ylabel('PLV')
    xlabel('Signal length [s]')
    xtickangle(45)
    grid on
    ylim([0 1])
    legend([theoFH intPerFH], 'Theoretical prediction', ['Integer num. ' ...
                        'of period'])
    set(gca, 'fontsize',fscp)

    box off
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength = [0.0155 0.035];

    hold on
    r2.x = -1.4;
    r1.y = 3.7;

    % putting subplot labbel
    text(r2.x, 1.25, 'D', 'HorizontalAlignment', 'center', 'FontWeight','bold', 'fontsize',vc.sblfz-2)
    text(r2.x, r1.y, 'A', 'HorizontalAlignment', 'center', 'FontWeight','bold', 'fontsize',vc.sblfz-2)
    text(7, r1.y, 'B', 'HorizontalAlignment', 'center', 'FontWeight','bold', 'fontsize',vc.sblfz-2)
    text(14.25, r1.y, 'C', 'HorizontalAlignment', 'center', 'FontWeight','bold', 'fontsize',vc.sblfz-2)

    %% the previous analysis, but demonstration for a 3 gamma_T
    gamma_Ts = [.75 .5 1]
    for kExample = 1 : 3
        subplot2(nR,nC, 1, kExample)
        gamma_T = gamma_Ts(kExample);
        hPSandOscParams.gamma_T = gamma_T;
        for k = 1 : hPSandOscParams.nSim
            tmpPLV(k) = cmpt_PLV_wPoissonSpkAndLinPhase(hPSandOscParams);
        end

        % plotting PLVs in a complex plane
        polarplot(angle(tmpPLV), abs(tmpPLV), '.', 'color',vc.sdc)
        hold on
        polarplot(angle(mean(tmpPLV)), abs(mean(tmpPLV)), 'ro', 'linewidth', mlw, 'markersize',ms)

        abs(mean(tmpPLV))
        % tgt : temporary variable for temporary ground truth PLV
        tgt = (2*pi*gamma_T*1i)^-1 * (exp(2*pi*gamma_T*1i) - 1);
        polarplot(angle(tgt), abs(tgt), ...
                  '+', 'linewidth',mlw, 'markersize',ms, 'color',vc.gtc)

        % do some adjustment on the figures (only visualization)
        tr = get(gca, 'rlim');
        trn = find_closestEvenNum(tr(2));
        rlim([0 trn])
        bfypp('PLV')
    end

end

function [PLV, linearPhase, varargout] = cmpt_PLV_wPoissonSpkAndLinPhase(hPSandOscParams)
% this function compute PLV with homogeneous Poisson spike train and perfectly linear phase

    T = hPSandOscParams.gamma_T / hPSandOscParams.f;
    t = 0 : hPSandOscParams.dt : T;

    % generate an complex exponential 
    cmplxValOsc = exp(1i * 2 * pi * hPSandOscParams.f * t );

    linearPhase = linspace(0, 2*pi*T*hPSandOscParams.f, numel(cmplxValOsc));
    varargout{1} = cmplxValOsc;
    t = (1:numel(cmplxValOsc)) * hPSandOscParams.dt;

    
    nOcs = 1; % it's not necessary, but could be use in case needed
    nUnit = 1; % it's not necessary, but could be use in case needed
    nSim = 1; % it's not necessary, but could be use in case needed
              % rateFact = 1;


    PLV = zeros(nOcs,nUnit,nSim);
    nSpk = zeros(nUnit,nSim);
    
    for ksim = 1:nSim
        for ktrial = 1:hPSandOscParams.nTrial
            spkTimes = ...
                rand(length(t),nUnit)<(ones(length(t),1)* hPSandOscParams.rate*hPSandOscParams.dt);
            % lfp = y;
            PLV(:,:,ksim) = PLV(:,:,ksim) + exp(1i*linearPhase) * spkTimes;
            nSpk(:,ksim) = nSpk(:,ksim)+sum(spkTimes,1)';
        end
    end
    
    PLV = PLV / nSpk;

end


% the following functions are only used for visualization purposes
function bfypp(rlabel)
% BeautiFY the Polar Plots

    thetaticks(0:45:315)

    rRange = get(gca, 'rlim')
    tn = rRange(2);

    rticks([0 tn/2 tn]);
    rlim([0 tn])
    ax = gca;

    rruler = ax.RAxis;

    ax.RAxisLocation = 22.5;

    ax.GridAlpha = .2;
    ax.LineWidth = 1.5;

    rruler = ax.RAxis;
    rruler.Label.String = rlabel;
end

function num = find_closestEvenNum(numIn)
% this function find the closet even number (to be used for r-ticks)
    tn = round(numIn * 100);
    if mod(tn, 2) == 0, 
        tn = tn; 
    else
        tn = tn + 1;
    end
    
    num = tn / 100;
end

