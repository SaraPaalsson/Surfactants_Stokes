function cb = colorlinesplot(x,y,c,figindex,specline,nwidth,speccolor,cmin,cmax,colorb)
% Colour line over x,y and z according to c.
%   - 2D plot specified by z=[]
%   - Plot in figure figindex, default 99.
%   - Linewidth specified by specline, default 1.
%   - Colormap specified by speccolor, default viridis.
% Updated 2016-10-07.
% Updated 2018-05-11 to add surf plot and return cb

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set parameters
if isempty(figindex)
    figindex = 99;
end

if isempty(specline)
    specline = 1;
end

if isempty(speccolor)
    viridisM=viridis();
    speccolor = viridisM;
end
if isempty(nwidth)
    nwidth = 0.025;
end

if isempty(colorb)
    colorb = 1;
end

if isempty(cmin)
    cmin = min(c);
end
if isempty(cmax)
    cmax = max(c);
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
my_sfigure(figindex);
cmap = colormap(speccolor);


% % change c into an index into the colormap
% % min(c) -> 1, max(c) -> number of colors
% cold = c;
% if cmax == cmin
%     c = ones(size(c));
% else
%     c = round(1+(size(cmap,1)-1)*(c - cmin)/(cmax-cmin));
% end

% Create and plot z interface:
zk = x + 1i*y;
zk2 = [zk; zk(1)];

if max(c) > 0    
    zkp = fft_diff(zk,[1 length(zk)]);
    n = -1i*zkp./abs(zkp);
    n = [n; n(1)];
    zkn = [zk2 zk2+nwidth*n]; %Create band inwards for surfactant concentration
    Sk = [c; c(1)];
    
    hs=surf(real(zkn),imag(zkn),zeros(size(zk2,1),2),[Sk Sk],'EdgeColor','interp') ;
    view(2) %Set view in X-Y plane
    hold on
end

plot(zk2,'k-','LineWidth',specline)
hold on

if colorb
    cb = colorbar;
    cmin = 0.9*cmin;
    cmax = 1.1*cmax;
    caxis([cmin cmax]);
else
    cb = 0;
end

end