% This example reproduces results from section 4.9 of the book. It
% simulates a series of dynamic charge-depleting UDDS profiles at 25 degC,
% starting at 95% SOC with output blending and model blending.
% 
% Copyright (©) 2024 The Regents of the University of Colorado, a body
% corporate. Created by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCCS). This work is licensed
% under a Creative Commons "Attribution-ShareAlike 4.0 International" Intl.
% License. https://creativecommons.org/licenses/by-sa/4.0/ 
% This code is provided as a supplement to: Gregory L. Plett and M. Scott
% Trimboli, "Battery Management Systems, Volume III, Physics-Based
% Methods," Artech House, 2024. It is provided "as is", without express or
% implied warranty. Attribution should be given by citing: Gregory L. Plett
% and M. Scott Trimboli, Battery Management Systems, Volume III:
% Physics-Based Methods, Artech House, 2024.         

clearvars; close all; clc;
if(~isdeployed),cd(fileparts(which(mfilename)));end
addpath(genpath('../'));

configXRA = 'defaultHRA.xlsx';   % xRA tuning values
inputFile = 'inputDynamic.xlsx'; % profile of current/temp vs time

% The following files correspond to the book notes:
cellParam   = 'cellNMC30.xlsx'; % parameter file name
ROMfile     = 'ROM/ROM_NMC30_HRA.mat'; % generated ROM
ROMoutFile1 = 'ROM/OUT_NMC30_DYNAMIC_MDLB.mat'; % ROM simulated output, mdlBlend
ROMoutFile2 = 'ROM/OUT_NMC30_DYNAMIC_OUTB.mat'; % ROM simulated output, outBlend
FOMfile     = 'FOM/FOM_NMC30_DYNAMIC.mph'; % generated FOM
FOMoutFile  = 'FOM/OUT_NMC30_DYNAMIC.mat'; % FOM simulated output

%% Generate ROM
% Note that ROM generation takes about one minute per SOC setpoint but
% happens one time only (most time is spent in C2D functionality)
try % Load cell parameters to "cellData" structure
  cellData = loadCellParams(cellParam);
catch
  error('Cell parameter file does not exist!');
end
  
if ~exist(ROMfile,'file') % then, need to create the ROM
  try % Load xRA tuning parameters to "xraData" structure
    xraData = loadXRA(configXRA);
    % uncomment the next line if you wish to see debug plots comparing the
    % ideal and ROM frequency responses
    % xraData.debug = 1;     
  catch
    error('xRA configuration file does not exist!');
  end
  
  ROM = genROM(cellData,xraData,'HRA');
  save(ROMfile,'ROM');
else % load the already-created ROM
  load(ROMfile,'ROM');
end

%% Simulate ROM
try % Load input current profile to "simData" structure
  simData = loadInput(inputFile);
catch
  error('Input current profile does not exist!');
end

% First, simulate using model blending
ROMout = simROM(ROM,simData,'mdlBlend'); ROMoutMDL = ROMout;
save(ROMoutFile1,'ROMout');
% Next, simulate using output blending
ROMout = simROM(ROM,simData,'outBlend'); ROMoutOUT = ROMout;
save(ROMoutFile2,'ROMout');

%% Generate and simulate FOM
if ~exist(FOMoutFile,'file') % then, need to create the ROM
  try % see if COMSOL LiveLink for MATLAB exists and is running
    import com.comsol.model.util.*
    % Generate FOM
    fprintf('Generating FOM...\n');
    FOM = genFOM(cellData);

    % Simulate FOM
    [FOM,FOMout] = simFOM(FOM,simData);
    mphsave(FOM,FOMfile);
    save(FOMoutFile,'FOMout');

    % uncomment next line to load FOM into COMSOL GUI
    % mphlaunch(mphtags(FOM));
  catch
    fprintf('COMSOL LiveLink for MATLAB does not appear to be running.\n');
    fprintf('Using precomputed data instead.\n');
    load(FOMoutFile,'FOMout');
  end
else
  load(FOMoutFile,'FOMout');
end

%% Plot results
% First, plot the model-blending results
plotResults(ROMoutMDL,FOMout,1,10)
% Next, plot the output-blending results
plotResults(ROMoutOUT,FOMout,1,10)
