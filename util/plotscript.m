% Script for plotting drops and their surfactant concentrations
% 
% Created 2018-05-11

close all; clear all; clc

defaultplotscript

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% If needed set which file to plot from
fileName = '../output/data/swissroll3F_outdata_t0p1.mat';
% fileName = '../output/data/swissroll3F_outdata_t2p6.mat';
% fileName = '../output/data/swissroll3F_outdata_t5p1.mat';
% fileName = '../output/data/swissroll3F_outdata_t10p1.mat';
% fileName = '../output/data/swissroll3F_outdata_t20p1.mat';
% fileName = '../output/data/swissroll3F_outdata_t30p1.mat';
% fileName = '../output/data/swissroll3F_outdata_t50p1.mat';
load(fileName)
z = simdata.drops.z;
idx = simdata.drops.idx;
rho = simdata.surfact.rho;

%% Plot interfaces with surfactant concentration on them

figure(1); clf
minc = inf; maxc = 0; 
minx = inf; maxx = 0;
miny = inf; maxy = 0;
for k=1:size(idx,1)
    zk = z(idx(k,1):idx(k,2)); 
    rhok = rho(idx(k,1):idx(k,2));

    %Plot interface with surfactant concentration in colour
    cb = colorlinesplot(real(zk),imag(zk),rhok,1,[],[],[],min(rhok),max(rhok),0);
    
    minc = min(minc,min(rhok));
    maxc = max(maxc,max(rhok));
    minx = min(minx,min(real(zk)));
    maxx = max(maxx,max(real(zk)));
    miny = min(miny,min(imag(zk)));
    maxy = max(maxy,max(imag(zk)));
end


% ylabel(cb,'Surfactant concentration','FontSize',20,'interpreter','latex')
% caxis([minc maxc])

if minx < 0
    minx = 1.2*minx;
else 
    minx = 0.8*minx;
end
maxx = 1.2*maxx;
if miny < 0
    miny = 1.2*miny;
else 
    miny = 0.8*miny;
end
maxy = 1.2*maxy;

xlim([minx maxx])
ylim([miny maxy]);

axis equal

xlabel('$x$'); 
ylabel('$y$');

% Paper swissroll
minc = 0; maxc = 6.5;
caxis([minc maxc])
axis([-1.9 1.1 -1.2 1.75])
% set(cb,'Visible','Off');
% set(cb,'Location','North')
set(gca,'Visible','Off')
grid off
box on

%% Plot x and y coordinates as function of alpha for drop k

k = 1; 
zk = z(idx(k,1):idx(k,2)); zk = [zk; zk(1)];
alpha = linspace(0,2*pi,length(zk))';

figure(2); clf
plot(alpha,real(zk),'k-','DisplayName','$x(\alpha)$'); 
hold on
plot(alpha,imag(zk),'k-.','DisplayName','$y(\alpha)$')
hl = legend('toggle');
set(hl,'FontSize',15,'interpreter','latex')
xlabel('$\alpha$')

%% Plot interface of drop k and parametrisation orientation

k = 1; 
zk = z(idx(k,1):idx(k,2)); 
zpk = fft_diff(zk,[1 length(zk)]);
zk = [zk; zk(1)];

t = zpk(1); 
n = -1i*t;

zn = z(1);

figure(3); clf
plot(zk,'k-');
hold on
quiver(real(zn),imag(zn),real(t)/abs(t),imag(t)/abs(t),0.2,'MaxHeadSize',0.8);
quiver(real(zn),imag(zn),real(n)/abs(n),imag(n)/abs(n),0.2,'MaxHeadSize',0.8);


%% Plot surfactant concentration vs alpha with the graph coloured according to rho

k = 1; 
rhok = rho(idx(k,1):idx(k,2));
rhok2 = [rhok; rhok(1)];
alpha = linspace(0,2*pi,length(rhok)+1)';

cmap = colormap('jet');

rhokp = fft_diff(rhok,[1 length(rhok)]);
n = -1i*rhokp./abs(rhokp);
n = [n; n(1)];

nwidth = 0.01;
    
rhokn = [rhok2-nwidth./rhok2 rhok2+nwidth./rhok2]; %Create band inwards for surfactant concentration
   
figure(4); clf
hs=surf([alpha alpha],rhokn,zeros(size(rhok2,1),2),[rhok2 rhok2],'EdgeColor','interp') ;
view(2) %Set view in X-Y plane
hold on




