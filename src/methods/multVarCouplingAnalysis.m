function [svdOut, varargout] = multVarCouplingAnalysis(spikeTrains_raw, lfpPhases_raw, varargin)
% this function do the multi-variate phase locking analysis and return the singular values
    

    %% prepare the spike train and LFPs
    [spikeTrains_allTrLong, lfpPhases_allTrLong, n, selectedUnits, unwhitenOpr] = ...
        prep_SpkLfpData(spikeTrains_raw, lfpPhases_raw, varargin{:});

    %% compute the coupling matrix
    [M] = cmpt_couplingMatrix(...
        full(spikeTrains_allTrLong), lfpPhases_allTrLong);

    %% apply SVD on coupling matrix
    [singularLfpVecs_raw, singularSpkVecs_raw, singularValues] = fctrz_couplingMatrix(M);
    
    svdOut.singularValues = singularValues;
    svdOut.couplingMatrix = M;
end


function [spikeTrains_allTrLong, lfpPhases_allTrLong, n, selectedUnits, ...
          unwhitenOpr] = prep_SpkLfpData(spikeTrains_raw, lfpPhases_input, varargin)

    %%
    n.Tr        = size(lfpPhases_input, 3);        % number of trials

    %% whitening the LFP (the continious signal)
    [tmpSig, witenOpr, unwhitenOpr] = ...
        whitenSignal(reshape(lfpPhases_input, size(lfpPhases_input, 1), []), NaN);
    % when is NaN, no rank reduction will be applied 
    
    lfpPhases_raw = reshape(tmpSig, size(tmpSig, 1), [], n.Tr);

    %% signal info
    n.LfpCh     = size(lfpPhases_raw, 1);        % number of LFP channels
    n.Sample    = size(lfpPhases_raw, 2);        % number of samples
    n.SpkUnit   = size(spikeTrains_raw{1}, 1);   % number of spiking units

    % all sample are selected
    selectedSamples = 1 : n.Sample;

    
    %% Prepare the data
    % no prep for LFP 
    lfpPhases = lfpPhases_raw;

    % preparing spike matrix
    for iTr = 1 : n.Tr
        spikeTrains{iTr} = spikeTrains_raw{iTr}(:, selectedSamples);
    end
    selectedUnits = [];

    %% Creat long spike trains and LFP 
    spikeTrains_allTrLong = cell2mat(spikeTrains);

    updated_nSample = numel(selectedSamples);
    lfpPhases_allTrLong = nan(n.LfpCh, n.Tr * updated_nSample);
    for iTr = 1 : n.Tr
        tmpRange = (iTr-1)*updated_nSample +1 : (iTr-1)*updated_nSample + updated_nSample;
        lfpPhases_allTrLong(:,tmpRange) = lfpPhases(:,:, iTr);
    end
end

function [whitenDataMat, varargout] = whitenSignal(dataMat_raw,proportion)
% select PCs accounting for proportion*100 pct of the variance and whiten
% the signal by projecting on those (this also reduced dimension)
% if NaN is passed, no reduction is applied
% inputs:
% * dataMat_raw: data matrix (channels/features x samples)
% * proportion: proportion <1 of variance to account for, defaut .99, or
% number of components to keep (integer>=1)
    if nargin<2
        proportion = .99;
    end
    N = size(dataMat_raw, 2);
    dataMat = bsxfun(@minus, dataMat_raw, mean(dataMat_raw, 2));
    covMat = dataMat * dataMat' / N;

    [U, D] = eig(covMat);
    [d,ind] = sort(diag(D),'descend');
    cumTrace = cumsum(d);
    if proportion<1
        n = find(cumTrace>=proportion*cumTrace(end));
        n = n(1);
    elseif isnan(proportion)
        n = size(diag(D), 1);
    else
        n = proportion; 
    end
    ind = ind(1:n);
    Ds = D(ind,ind);
    Us = U(:,ind);

    whitenDataMat = diag(diag(Ds) .^ -.5) * Us' * dataMat;

    whitenOpr = diag(diag(Ds) .^ -.5) * Us';
    whitenInvOpr = Us * diag(diag(Ds) .^ .5);
    varargout{1} = whitenOpr;
    varargout{2} = whitenInvOpr;
    varargout{3} = mean(dataMat_raw, 2);
end

function [couplingMatrix varargout] = cmpt_couplingMatrix(spkTrain, lfpPhase, varargin)
% compute the coupling matrix

    %%
    nCh = size(lfpPhase, 1);

    %%
    % *** we might need to run some check here:
    % the phases are in range ([-pi pi])
    % user want to use the envelope

    if isreal(lfpPhase)
        F = exp(1i * lfpPhase);
    else
        F = abs(lfpPhase) .* exp(1i * angle(lfpPhase));
    end

    varargout{1} = sum(spkTrain, 2);


    normalizingFactor = sum(spkTrain,2 ).^ -.5;

    % take care of the case with zero number of spikes
    normalizingFactor(isinf(normalizingFactor))= nan;
    normalizingMatrix = repmat(normalizingFactor', nCh, 1);
    
    tmp_couplingMatrix = F * spkTrain';    
    couplingMatrix = ( normalizingMatrix .* abs(tmp_couplingMatrix) ) .* exp(1i * angle(tmp_couplingMatrix));


end

function [singularLfpVecs, singularSpkVecs, singularValues] = fctrz_couplingMatrix(M)
% [singularLfpVecs, singularSpkVecs, singularValues] = fctrz_couplingMatrix(M)
% apply SVD
    
    %% SVD decompostion
    if sum(isnan(M(:))) > 0 % means you had some units without any spike
        warning('You had units with no spikes (at least in some trials). These units were excluded from the analysis')
	[M_cleaned, ~, valIndex] = remove_NaNinMat(M, 'fullRowOrCol', 'col');
        [singularLfpVecs, D, V_cleaned] = svd(M_cleaned); 
        singularSpkVecs = nan(size(M, 2));
        singularSpkVecs(valIndex.c,1:numel(valIndex.c)) = V_cleaned;
    else % normal procedure
        [singularLfpVecs, D, singularSpkVecs] = svd(M); 
    end

    % if M is one-dimentional the diag will turn in into a matrix
    if any(size(M) == 1)
        singularValues = D;
    else
        singularValues = diag(D);
    end
end

function [mat_nanRemoved varargout] = remove_NaNinMat(mat, varargin)
% [mat_nanRemoved nans vals] = remove_NaNinMat(mat, NaNdistributionType, removeType)
% remove NaN values from a matrix
%
%
% EXAMPLE:
%
%
% ------
% Input:
% 1     mat: matrix to be cleaned
%
% Output:
% 1     mat_nanRemoved: clean matrix
%
% ------
% see also isnan
% ------
% potential improvments:
% (1) extending beyond 2D matrix
% (2) maybe write a new function which to the same thing with inf, empty
% guys collectivly
% ------
% Code Info:
%   creation: 2015-06-25 by SS (codes@shervinsafavi.org)
%   modification:
%       $ 201?-??-?? ?

    %% handle optional inputs (varargin): NaNdistributionType, removeType
    optionalVariables.NaNdistributionType = []; optionalVariables.removeType = [];
    defaultValues{1} = 'sparse'; defaultValues{2} = 'both';
    optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

    [nR nC] = size(mat);

    switch optionalVariables.NaNdistributionType
      case 'sparse'
        
        [tmp_nans.r tmp_nans.c] = find(isnan(mat));
        [tmp_vals.r tmp_vals.c] = find(~isnan(mat));
        
        nans.r = unique(tmp_nans.r);
        nans.c = unique(tmp_nans.c);
        
        vals.r = setdiff((1:nR), nans.r);
        vals.c = setdiff((1:nC), nans.c);

        switch optionalVariables.removeType
            
          case 'row'
            mat_nanRemoved = mat(vals.r, :);
          case 'col'
            mat_nanRemoved = mat(:, vals.c);
          case 'both'
            mat_nanRemoved = mat(vals.r, vals.c);
          otherwise
        end
        
      case 'fullRowOrCol'
        
        logicalIndex_mat = isnan(mat);
        
        % check which rows are fully NaN
        nNaNinRows = sum(logicalIndex_mat, 2);
        nans.r = find(nNaNinRows == nC);
        vals.r = setdiff((1:nR), nans.r);
        
        % check which columns are fully NaN
        nNaNinCols = sum(logicalIndex_mat, 1);
        nans.c = find(nNaNinCols == nR);
        vals.c = setdiff((1:nC), nans.c);
        
        switch optionalVariables.removeType
          case 'row'
            mat_nanRemoved = mat(vals.r, :);
          case 'col'
            mat_nanRemoved = mat(:, vals.c);
          case 'both'
            nans.bothRowAndCol = union(nans.r, nans.c);

            vals.r = setdiff((1:nR), nans.bothRowAndCol); 
            vals.c = setdiff((1:nC), nans.bothRowAndCol);    
            mat_nanRemoved = mat(vals.r, vals.c);
          
          otherwise
            
        end
      otherwise
        
    end
    varargout{1} = nans;
    varargout{2} = vals;
end