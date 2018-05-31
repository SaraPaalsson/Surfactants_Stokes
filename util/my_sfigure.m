function h = my_sfigure(varargin)
% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%
% Updated for figure title and size, Sara Paalsson 2016.
%
% See also figure

if nargin>=1 
    h = varargin{1};
	if ishandle(h)
		set(0, 'CurrentFigure', h);
    else
        switch nargin
            case 1
                h = figure(h);
            case 3
                h = figure(varargin{2},varargin{3});
            case 5
                h = figure(varargin{2},varargin{3},varargin{4},varargin{5});
        end
	end
else
	h = figure;
end
