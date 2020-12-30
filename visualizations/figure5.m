function varargout = figure5()
% figure4()
% this function plot the figure 5 of following paper:
% [From univariate to multivariate coupling between continuous signals and point processes: a mathematical framework, S.Safavi, N. K. Logothetis and M. Besserve. ArXiv 2020](https://arxiv.org/abs/2005.04034)

    clf

    %% assign parameters

    % fix the seed of the random number generator to get consistent figure
    rng(2);

    % assign visualization parameters 
    vc = get_vizConventions();
    adjustFigAppearnce(vc.f5);
    lw = 2.5;
    
    % assign parameters for generating homogeneous Poisson spikes and accompanying oscillation 
    % hPSandOscParams: homogeneous Poisson spikes and [accompanying] oscillation parameters 

    % nRel = 100;
    nRel = 2;    

    % case without coupling 
    caseName = 'noCoupling';
    couplingStrength = 0;
    vmCMwMPparams = assign_vmCMwMPparams(nRel, couplingStrength);
    % [svdOut_noCoupling] = wrapper_multiVarCouplingSimulation(vmCMwMPparams, caseName);
    [tmpSvdOut, vmCMwMPparams] = wrapper_multiVarCouplingSimulation(vmCMwMPparams, caseName);
    for kun = 1 : numel(vmCMwMPparams.unitNums)
        svdOut(kun).(caseName) = tmpSvdOut(kun).(caseName);
    end
    
    % case with coupling 
    caseName = 'withCoupling';
    couplingStrength = .15;
    vmCMwMPparams = assign_vmCMwMPparams(nRel, couplingStrength);
    [tmpSvdOut, vmCMwMPparams] = wrapper_multiVarCouplingSimulation(vmCMwMPparams, caseName);
    for kun = 1 : numel(vmCMwMPparams.unitNums)
        svdOut(kun).(caseName) = tmpSvdOut(kun).(caseName);
    end

    viz_EmpTheoDistributions(svdOut, vmCMwMPparams, vc)    
    
    % optional output is the 'svdOut' of both cases
    varargout{1} = svdOut;

end

function vmCMwMPparams = assign_vmCMwMPparams(nRel, couplingStrength)
% this function assign the parameters for the von Mises Coupling Model with Multiple Populations 

    %% different choices for the number of point processes (spikes) 
    vmCMwMPparams.unitNums = [10 50 90];
    nUnitNum = numel(vmCMwMPparams.unitNums);
    
    %% assigning signal parameters 
    for iun = 1 : nUnitNum % loop on different choices of number of units
        vmCMwMPparams.signalParams = struct ...
            ( ...
                'nCh',      NaN, ...         % number of LFP channel (oscillatory signals)
                'nUnit',    NaN, ...         % number of spiking units
                'SF',       1e3, ... % Hz    % sampling frequency
                'nTr',      10, ...          % number of trials
                'signalLength', ...  % sec   % duration of the signals
                11 ...
                );
    end

    
    % frequencies of oscillatory components  
    vmCMwMPparams.allOscFreq = [11  12 13 14 15];
    % baseline mixing coefficient
    vmCMwMPparams.mixingBaseLineCoef = 0;
    % number of realization of the simulation
    vmCMwMPparams.nRel = nRel;

    % choice for the number of oscillatory signals (LFPs)
    vmCMwMPparams.signalParams.nCh = 100;

    %% firing rate or intensity function 
    % this value is used for all spiking populations 
    vmCMwMPparams.aveFR = 20;

    vmCMwMPparams.couplingStrength = couplingStrength;
    
    % assigning the parameters of global dynamics
    clear globalDynamicsParams

    for iun = 1 : nUnitNum % loop on different choices of number of units
        vmCMwMPparams.globalDynamicsParams(iun) = struct ...
            (...
                'oscFreq',  vmCMwMPparams.allOscFreq, ...  % Hz
                'lfpPhaseNoise_kappa', 10, ...   % concentration parameter of noise on LFP phase
                'whiteNoise_sigma', ...      % variance for white noise (not used)
                0 ...  
                );
    end

    %% parameters for dividing to multiple populations
    % noc : number of clusters - used for different frequencies of oscillation 
    for iun = 1 : nUnitNum
        noc = numel(vmCMwMPparams.globalDynamicsParams(iun).oscFreq);
        vmCMwMPparams.globalDynamicsParams(iun).nFreqComp = numel(vmCMwMPparams.globalDynamicsParams(iun).oscFreq);
        vmCMwMPparams.globalDynamicsParams(iun).oscComps = ones(1, noc);
    end

    %% filter used for LFP preprocessing
    vmCMwMPparams.filterOrder = 2;
    % a form of bias correction on phase that is not used in the analysis
    vmCMwMPparams.nIteCent = 0; 
    vmCMwMPparams.halfFilterWidth = 5;
    iFreq = 1; % frequency of LFP oscillation is used for spike-LFP coupling analysis
    % adjusting the centeral frequency for filtering
    vmCMwMPparams.freqCenter = vmCMwMPparams.globalDynamicsParams(iun).oscFreq(iFreq);

end

function viz_EmpTheoDistributions(svdOut, vmCMwMPparams, vc)
% this function visualize the empirical and theoretical distribution of eigenvlaues
    
    nR = vc.f5.nR;
    nC = vc.f5.nC;

    fn = fieldnames(svdOut(1));
    nCase = numel(fn);
    
    % loop on simulation with and without coupling 
    for kcase = 1 : nCase
        columnID = kcase;
        caseName = fn{kcase};
        
        if kcase == nCase, 
            legendFlag = 1;
        else 
            legendFlag = 0;
        end
        
        nUnitNum = numel(svdOut);
        
        for iun = 1 : nUnitNum
            subplot2(nR,nC, iun, columnID);

            %%
            clear msd td_all
            tmp2ods = [];
            
            % loop on different realization of the simulations
            for iRel = 1 : vmCMwMPparams.nRel
                % converting the singular values to eigenvalues 
                tmp = ...
                    (svdOut(iun).(caseName){1, iRel}.singularValues) .^ 2 ...
                    / vmCMwMPparams.signalParams(iun).nUnit; 
                tmp2ods = [tmp2ods tmp(:)];
            end
            tods = tmp2ods(:);
            
            % histogram bins for later
            histAx = linspace(0, max(tods)+1, 50);

            %% plotting the theoreitical and empirical distributions 
            td = hist(tods, histAx) / numel(tods);
            b = bar(histAx, td, 'facecolor',.5*ones(1,3), 'EdgeColor',.4*ones(1,3), ...
                    'BarWidth',.9)

            % aspect ratio of the coupling matrix
            c = ...
                size(svdOut(iun).(caseName){1, iRel}.couplingMatrix, 1) ...
                / size(svdOut(iun).(caseName){1, iRel}.couplingMatrix, 2);
            lambda = (1 + c^.5) ^ 2;

            % variance of individual entries of coupling matrix is 1 given normalizations implemented in earlier
            s = 1; % variance

            % theoretical left and right boundaries of the spectra distribution 
            a=(s^2)*(1-sqrt(c))^2;
            b=(s^2)*(1+sqrt(c))^2;

            % grid for singular values to find the empirical distribution 
            svGrid = 0 : diff(histAx(1:2)) : lambda + diff(histAx(1:2));
            
            % theoretical PDF
            ft=@(svGrid, a,b,c) (1./(2*pi* svGrid * c * s^(2))).*sqrt((b - svGrid).*(svGrid - a));
            F = real(ft(svGrid,a,b,c)); 
            % for the values outside domain defined real value will be zero
            pdf = F/(sum(F));

            hold all
            plot(svGrid, pdf, 'color',vc.gtc, 'LineWidth', vc.f5.lw);
            if iun == 1    
                xline(lambda, 'b', 'linewidth',2);
            else
                xline(lambda, 'b', 'linewidth',2);
            end
            axis tight; box off
            ax = gca;
            set(ax, 'TickLength', [0.04 0.035])
            set(ax, 'TickDir', 'out')

            %%

            if iun == nUnitNum
                xlabel('Eigenvalue')
            end

            axes('Position', vc.f5.axisLoc{columnID}(iun, :));
            box on; 
            b = bar(histAx, td, 'facecolor',.5*ones(1,3), 'EdgeColor',.4*ones(1,3), ...
                    'BarWidth',.9)
            hold all
            plot(svGrid, pdf, 'color',vc.gtc, 'LineWidth', vc.f5.lw);
            axis tight
            box on
            ax = gca;
            set(ax, 'TickLength', [0.04 0.035]);
            set(ax, 'ytick', []);

            ylim([0 .02]);
            xlim(vc.f5.xlimSel(iun, :));
            
        end
        
        % *** set the same lim for the other subplot

        if legendFlag == 1
            h = legend(...
                'Empirical distribution', 'Theoretical prediction', ...
                'Theoretical significance threshold', 'location','north');

            pos = get(h, 'Position');
            set(h, 'box', 'off');
            set(h, 'fontsize', 9);
            posx = .59; 
            posy = .92;
            set(h,'Position',[posx posy pos(3) pos(4)]);

        end

    end
end


