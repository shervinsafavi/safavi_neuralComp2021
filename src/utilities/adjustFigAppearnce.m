function adjustFigAppearnce(vizParams, varargin)

if nargin < 3
    fh = gcf;
end

fh.Color        = [1 1 1]
fh.Units        = 'centimeters';
fh.Position(3)  = vizParams.w;
fh.Position(4)  = vizParams.h;
