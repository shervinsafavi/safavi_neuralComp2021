function vc = get_vizConventions()
% get_vizConventions()
% this function give the parameters are needed only in visualization (visualization conventions)

%% general settings
vc.gtc = [0 1 0];
vc.sdc = .4 * ones(1,3);
vc.spkThings = 'r';
vc.lfpThings = 'b';

% size subplot labels
vc.sblfz = 15;

% for legend font size
vc.lfs = 14;


%% figure 1
vc.f1.w = 16.3072;
vc.f1.h = 8.2176;
vc.f1.mc = [.5 0 .5]

%% figure 2
vc.f2.w = 19.4235;
vc.f2.h = 24.4608;
vc.f2.nR = 1+2+1+3+1;
vc.f2.nC = 2+2;
% hfi : half figure index
vc.f2.hfi = (vc.f2.nC) / 2;
vc.f2.tlw = 2.8;         % trace line width
% specify the location of an extra axis 
tyLoc = .05
vc.f2.r4pos = [ ...
    0.13 tyLoc 0.156648936170213 0.0760103626943005; ...
    0.336117021276596 tyLoc 0.156648936170213 0.0760103626943005; ...
    0.542234042553192 tyLoc 0.156648936170213 0.0760103626943005; ...
    0.748351063829787 tyLoc 0.156648936170213 0.0760103626943005];
vc.f2.r3pos = [0.13 0.240569948186528 0.775 0.287150259067357];
vc.f2.fsi = 10;          % font size inset
vc.f2.ncc = [1 .65 0];   % no coupling color
vc.f2.wcc = [.7 0 .7];   % with coupling color
vc.f2.ll = 0.682;        % legend loc
vc.f2.imts = 0.05;       % 
vc.f2.fsm = 12;          % font size main

%% figure 3
vc.f3.w = 16.3072;
vc.f3.h = 8.2176;
vc.f3.nR = 2;
vc.f3.nC = 3;
