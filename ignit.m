tmpPath = which('ignit');
[phd, ~, ~] = fileparts(tmpPath); % phd : project home directory
                                  
% add necessary routines
addpath(fullfile(phd, 'src'));

