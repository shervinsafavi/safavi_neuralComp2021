tmpPath = which('ignit');
[phd, ~, ~] = fileparts(tmpPath); % phd : project home directory
                                  
% add necessary routines
addpath(genpath(fullfile(phd, 'src')));

clear phd tmpPath