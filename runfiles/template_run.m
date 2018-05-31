% Template run-file
% Updated 2018-03-07

close all;
clear;
% clc

for j=1:2
    try startup_script
        % We can run the script and should be were we should be. Awesome.
        disp('Startup_script complete.')
        break
    catch
        if j==1
            disp('No startup_script found. Stepping up one step and trying again...');
            cd ..
        else
            disp('Cannot find startup_script!')
            return
        end
    end
end


disp('***************************************************************');
disp('Setting up inparameters...');

% Define indata file:
fileInV = {'indata/template_indata'};

% In case of restart, define other indata file:
% fileInV = {'indata/template_restart'};

% Define outdata file
fileOutV = {'template_outdata'};


for i=1:length(fileInV)
    fileIn = [fileInV{i} '.m'];
    fileOut = ['output/data/' fileOutV{i}];

    try load(fileOut)
        disp(['File ' fileOut ' successfully loaded'])
    catch
        disp(['File ' fileOut ' does not exist. Runs simulation.'])
        [final] = stokes_surf(fileIn,fileOut,1);
        fileOutFinal = ['output/final_' fileOutV{i} '_t' num2str(final.t)];
        fileOutFinal(fileOutFinal == '.') = '';
        save(fileOutFinal,'final')
    end
   
end
