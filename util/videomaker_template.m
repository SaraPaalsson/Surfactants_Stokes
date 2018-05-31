% This scrip makes a video from files in the outdata/data folder.
%
% Creation data: 2018-04-10
% Updated: 2018-05-11

close all; clear; clc


defaultplotscript
startup_script

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% File info

fileName = 'swissroll3F_outdata';   %Filename to make video of

tvec = (0.1:0.1:80.1)'; %Array containing save times for outdata

rhoInt = [0 6.5]; %Set max an min of surfactants for colourbar

txtC = [-1.8 -1.1];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% To create a movie or not:
ifVideo = 0;

vQuality = 80; %Videoquality:
vFormat = 'mp4';

vName = 'mymovie';


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Filepath to take save files from

filePath = ['output/data/' fileName];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set axis and colour of plot
myaxes = [-2 1.5 -1.3 1.7];
myC = [200 200 200]/256;
mSize = 15;
noff = 0.03;
myA_surf = 0.25;
myA_nosurf = 0.75;
lw = 1.5; %linewidth for black lines (i.e. for drops without surfactants

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Define colormap
cmap = colormap(viridis);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create video object
if ifVideo
    switch vFormat
        case 'mp4'
            videoName = [vName '.mp4'];
            v = VideoWriter(videoName,'MPEG-4');
    end
    v.Quality = vQuality;
    open(v);
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Make movie

% Go through time and capture frames
for j=1:length(tvec)
    
    fT = [filePath '_t' num2str(tvec(j))]; %Define file for time tvec(j)
    fT(fT == '.') = 'p'; %Remove decimal points
    load(fT) %Load file
    
    disp(['Plot for time ' num2str(tvec(j)-0.1) ])

    idx = simdata.drops.idx;
    
    my_sfigure(1);
    hold off
    for k = 1:(simdata.drops.nbrDrops)
        zk = simdata.drops.z(idx(k,1):idx(k,2));
        zkp = fft_diff(zk,[1 length(zk)]);
        zk = [zk; zk(1)];
        zkp = [zkp; zkp(1)];
        n = -1i*zkp./abs(zkp);
        
        zkn = [zk zk+noff*n]; %Create band inwards for surfactant concentration
        
        ff = fill(real(zk),imag(zk),myC);
        hold on
        
        if simdata.surfact.rhomass0(k) ~= 0
            % Drop k is surfactant covered
            
            set(ff,'FaceAlpha',myA_surf);
            
            Sk = simdata.surfact.rho(idx(k,1):idx(k,2));
            Sk = [Sk; Sk(1)];
            
            hs=surf(real(zkn),imag(zkn),zeros(size(zk,1),2),[Sk Sk],'EdgeColor','interp') ;
            view(2) %// view(0,90)          %// set view in X-Y plane
            plot(zk,'k-','LineWidth',lw)
            
            disp(['For drop ' num2str(k) ', rhomax = ' num2str(max(Sk))]);
            
        else
            set(ff,'FaceAlpha',myA_nosurf);
            
            plot(zk,'k-','LineWidth',lw)
        end
        
    end
    
    set(gca,'Visible','Off')
    sfigure(1);
    axis equal
    axis(myaxes)
    h = colorbar;
    ylabel(h,'Surfactant concentration','FontSize',15);
    caxis(rhoInt)
    txtstr = ['t = ' num2str(tvec(j)-0.1)];
    text(txtC(1),txtC(2),txtstr,'Fontsize',20)

    drawnow

    if ifVideo
        F = getframe(gcf);
        writeVideo(v,F);
    else
        pause(0.01)
    end

    % Compute smallest distance between drops
    Dminvec(j) = findClosestDist(simdata.drops.z,idx);
    Dminvec(j)
end

min(Dminvec)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Close video object
if ifVideo
    close(v)
end

disp('Done!')
