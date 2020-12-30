function figure2()
% figure2()
% this function plot the figure 2 of following paper:
% [From univariate to multivariate coupling between continuous signals and point processes: a mathematical framework, S.Safavi, N. K. Logothetis and M. Besserve. ArXiv 2020](https://arxiv.org/abs/2005.04034)

    clf

    %% assign parameters

    % fix the seed of the random number generator to get consistent figure
    rng(2);

    % assign visualization parameters 
    vc = get_vizConventions();
    adjustFigAppearnce(vc.f2);
    hfi = vc.f2.hfi;

    tmpidx = 1 : 3000; % indices of time points used for plotting
    nSpt   = 30;       % number of spike trains used for plotting
    rosId  = [0 2];    % row off-set ID, only uysed for visualization purposes


    % assign parameters von Mises coupling model

    vmCMparams.T      = 5;         % simulation length
    vmCMparams.dt     = .001;      % simulation step
    vmCMparams.t      = ...        % simulation time steps%
        0 : vmCMparams.dt : vmCMparams.T;      
    vmCMparams.rate   = 10;        % event rate
    vmCMparams.nTrial = 50;        % number of trials in each simulation
    vmCMparams.nSim   = 50;        % number of simulation 
                                   % *** last to needd to change to 5000 for the final run
    vmCMparams.f      = 1;         % oscillation frequency of continious signal
    
    couplingKappas    = [0 .5];     % coupling strengths of both models

    % in this simulation the frequency of oscillation is assumed to be 1 Hz

    nModel = numel(couplingKappas);

    
    %% simulate von Mises coupling model - no coupling case

    for kmodel = 1 : nModel
        
        ros = rosId(kmodel);
        
        % kappa = couplingKappa;
        vmCMparams.kappa = couplingKappas(kmodel);
        [PLV, allspt] = smlt_vonMisesCouplingMocel(vmCMparams);
        % storing both PLVs for later usage
        PLVs{kmodel} = PLV;

        %% ploting

        % raster plots
        subplot2d(vc.f2.nR,vc.f2.nC, -1+[2 3], ros+[1 hfi]);
        plot_eventRaster(full(allspt(1:nSpt, tmpidx)), {'color', vc.f2.ncc, 'linewidth',1.6})
        axis off; box off

        % add figure labels
        % put_f2_subplotsLabels(vc)

        % modulators
        subplot2d(vc.f2.nR,vc.f2.nC, 3, ros+[1 hfi]);
        plot(mean(allspt(:,tmpidx))/vmCMparams.dt, 'color', vc.sdc)
        hold on
        plot(vmCMparams.rate * exp(vmCMparams.kappa * cos(2 * pi * vmCMparams.t(tmpidx)')), ...
             'color',vc.f2.ncc, 'linewidth',vc.f2.tlw);
        axis off; box off

        % oscillatory signal 
        subplot2d(vc.f2.nR,vc.f2.nC, 4, ros+[1 hfi])
        tmpLfp = exp(1i * 2 * pi * vmCMparams.t);
        plot(real(tmpLfp(:,tmpidx)), 'color', vc.lfpThings, 'linewidth',vc.f2.tlw)
        ylim([-2.5 1])
        axis off
        box off

        % empirical and theoretical distribution 
        [~, tdfh] = viz_EmpTheoDists(PLV, vc, vmCMparams, ros);

        % add legened
        if kmodel == 2
            h = legend(tdfh, 'Theoretical prediction');
            pos = get(h,'Position')

            % top right
            posx = vc.f2.ll;
            posy = 0.15;
            % posy = 0.10;

            set(h,'Position',[posx posy pos(3) pos(4)]);
        end
    end

    
    %% plot of the empirical distribution in the complex plane
    % run v5_plotScatterHist.m
    viz_EmpTheoDists_2D(vc, PLVs, vmCMparams, couplingKappas)

end

function [fh, tdfh] = viz_EmpTheoDists(PLV, vc, vmCMparams, ros)

    %%
    distRange = .35
    histAx = linspace(-distRange, distRange, 15);

    %% real
    sPLV = ...
        sqrt(vmCMparams.nTrial) ...
        * (PLV - besseli(1, vmCMparams.kappa) / besseli(0,vmCMparams.kappa));
    distPLV = hist(real(sPLV), histAx);
    alpha_T = vmCMparams.T * vmCMparams.f;

    tmpSigma = ...
        ((besseli(0, vmCMparams.kappa) + besseli(2, vmCMparams.kappa)) ...
         / ((2 * vmCMparams.rate * vmCMparams.T * besseli(0, vmCMparams.kappa)^2))) ^ .5;

    fh.r4(ros+1) = axes('Parent',gcf, 'Position', vc.f2.r4pos(ros+1,:));
    % run v3_plotIndFit_core.m
    plotIndFit_core(vmCMparams, distPLV, histAx, vc, tmpSigma)
    xlabel('Real part')
    set(gca, 'fontsize', vc.f2.fsi)
    ax = gca;
    ax.YAxis.Exponent = 2;

    if ros == 0
        ylabel('Frequency');
    end

    %% imaginary
    sPLV = sqrt(vmCMparams.nTrial) * (PLV);
    distPLV = hist(imag(sPLV), histAx);
    alpha_T = vmCMparams.T * vmCMparams.f;
    tmpSigma = ...
        ((besseli(0, vmCMparams.kappa) - besseli(2, vmCMparams.kappa)) ...
         / ((2 * vmCMparams.rate * vmCMparams.T * besseli(0, vmCMparams.kappa)^2))) ^ .5;

    fh.r4(ros+2) = axes('Parent',gcf, 'Position', vc.f2.r4pos(ros+2,:));
    tdfh = plotIndFit_core(vmCMparams, distPLV, histAx, vc, tmpSigma);
    xlabel('Imaginary part')
    set(gca, 'fontsize', vc.f2.fsi)
    ax = gca;
    ax.YAxis.Exponent = 2;
end

function tdfh = plotIndFit_core(vmCMparams, distPLV, histAx, vc, tmpSigma)
    
    if vmCMparams.kappa == 0,  
        distRange = .01; 
        tc = vc.f2.ncc;
    else
        tc = vc.f2.wcc; 
    end

    bar(histAx, distPLV, 'FaceColor',tc, 'Facealpha',.9, 'edgecolor','none');
    hold on
    dAx = diff(histAx([1,2]));


    tdfh = plot(histAx, ...
                (1 / (tmpSigma * sqrt(2*pi))) * ...
                exp(-(histAx / tmpSigma).^2/2)* vmCMparams.nSim * dAx, ...
                'linewidth',2, 'color',vc.gtc);

    grid on
    mr = max(distPLV);
    ylim([0 mr*1.2]);

end

function viz_EmpTheoDists_2D(vc, PLVs, vmCMparams, couplingKappas)
    axes('Parent',gcf,...
         'Position', vc.f2.r3pos);


    %% assign PLVs for both models
    PLV_noCoupling = PLVs{1};
    PLV_wCoupling = PLVs{2};

    kappa = couplingKappas(1);
    plvgt = besseli(1,kappa) / besseli(0,kappa);
    ncdfh = plot(real(PLV_noCoupling), imag(PLV_noCoupling), '.', 'color', vc.f2.ncc)
    hold on
    gtfh = plot(plvgt, 0, '+', 'color',vc.gtc);

    kappa = couplingKappas(2);
    plvgt = besseli(1,kappa) / besseli(0,kappa);
    wcdfh = plot(real(PLV_wCoupling), imag(PLV_wCoupling), '.', 'color', vc.f2.wcc)
    hold on
    gtfh = plot(plvgt, 0, '+', 'color',vc.gtc);
    tyr = .1*[-.2 2]; % tmp y range
    ylim(tyr)
    yticks([0 max(tyr)/2 max(tyr)])
    xlim([-.02 .4])

    ylabel('Imaginary part')
    xlabel('Real part')

    set(gca, 'fontsize', vc.f2.fsm)

    box off
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength = [0.0155 0.035];

    %% legend
    h = legend([ncdfh wcdfh gtfh], ...
               'Simulation without coupling', 'Simulation with coupling', ...
               'Theoretical PLV', 'location','north')

    pos = get(h,'Position')


    % top right
    posx = vc.f2.ll - .122;
    posy = 0.485;

    set(h,'Position',[posx posy pos(3) pos(4)]);


    %% plot inset

    tyloc = .27; % temp y loc

    axes('Position',[.22 tyloc .2 .2])
    PLV = PLV_noCoupling;
    vmCMparams.kappa = couplingKappas(1);
    plotScatterHist_inset(PLV, vmCMparams, vc);
    axis square
    box on

    axes('Position',[.68 tyloc .2 .2])
    PLV = PLV_wCoupling;
    vmCMparams.kappa = couplingKappas(2);
    plotScatterHist_inset(PLV, vmCMparams, vc)

    axis square
    box on

end

function plotScatterHist_inset(PLV, vmCMparams, vc, couplingKappas)
    tv =  vmCMparams.T * vmCMparams.rate * sqrt(vmCMparams.nTrial) * ...
          (PLV - besseli(1, vmCMparams.kappa) / besseli(0, vmCMparams.kappa)) ...
          / sqrt(vmCMparams.T * vmCMparams.rate / 2) / 1;

    tv_wxms = PLV;

    plvgt = besseli(1, vmCMparams.kappa) / besseli(0,vmCMparams.kappa);

    if vmCMparams.kappa == 0,  
        distRange = .01; 
        tc = vc.f2.ncc;
    else
        tc = vc.f2.wcc; 
    end


    plot(real(PLV), imag(PLV), '.', 'color', tc)
    hold on
    plot(plvgt, 0, '+', 'color',vc.gtc, 'markersize', 10, 'linewidth',2);
    ax = gca;
    ax.XAxis.Exponent = -2;
    ax.YAxis.Exponent = -2;
    if vmCMparams.kappa > 0
        xtickformat('%,.2f')
    end

    set(ax, 'TickLength', [vc.f2.imts 0.035])

    set(gca, 'fontsize', vc.f2.fsi)

    ylabel('Imaginary part')
    xlabel('Real part')

    grid on

    %%
    if vmCMparams.kappa == 0,  
        distRange = max(abs((PLV))) * 1.1;
        distLimit = distRange * [-1 1];
        xlim(distLimit)
        ylim(distLimit)

    end


    if vmCMparams.kappa == .5, 
        distRange = max(abs((PLV))) * 1.1;
        distLimit = distRange * [-1 1];
        xlim(plvgt + distLimit)
        ylim(distLimit)
    end

    yticks([distLimit(1) 0 distLimit(2)]);
    xticks(plvgt + [distLimit(1) 0 distLimit(2)]);

end

function put_f2_subplotsLabels(vc)
    ros = 0;

    vc.f2.c1.x = -550;
    vc.f2.r1.y = -4;
    vc.f2.c2.x = 3200;
    vc.f2.r3.y = 127;

    vc.f2.sp11.sblfz = vc.sblfz;
    vc.f2.sp11.x = vc.f2.c1.x;
    vc.f2.sp11.y = vc.f2.r1.y;
    insert_spLabel(vc.f2.sp11, 'A')

    vc.f2.sp12.sblfz = vc.sblfz;
    vc.f2.sp12.x = vc.f2.c2.x;
    vc.f2.sp12.y = vc.f2.r1.y;
    insert_spLabel(vc.f2.sp12, 'B')


    vc.f2.sp21.sblfz = vc.sblfz;
    vc.f2.sp21.x = vc.f2.c1.x;
    vc.f2.sp21.y = 63;
    insert_spLabel(vc.f2.sp21, 'C')


    vc.f2.sp31.sblfz = vc.sblfz;
    vc.f2.sp31.x = vc.f2.c1.x;
    vc.f2.sp31.y = vc.f2.r3.y;
    insert_spLabel(vc.f2.sp31, 'D')


    vc.f2.sp32.sblfz = vc.sblfz;
    vc.f2.sp32.x = vc.f2.c2.x;
    vc.f2.sp32.y = vc.f2.r3.y;
    insert_spLabel(vc.f2.sp32, 'E')
end

function insert_spLabel(fvc, label)
% insert_spLabel(fvc, label)
    hold on
    text(fvc.x, fvc.y, label, 'HorizontalAlignment', 'center', 'FontWeight','bold', 'fontsize',fvc.sblfz)
end

