% This example reproduces results from section 4.10 of the book. It
% simulates the CC/CV and CP/CV charging-example results.
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

cellParam = 'cellNMC30.xlsx'; % parameter file name
configXRA = 'defaultHRA.xlsx';  % xRA tuning values
ROMfile   = 'ROM/ROM_NMC30_HRA.mat'; % generated ROM
ROM_CCCV  = 'ROM/ROMoutB_CCCV.mat'; % simulated ROM 
ROM_CPCV  = 'ROM/ROMoutB_CPCV.mat'; % simulated ROM 

%% Generate ROM if not already done
% Note that ROM generation takes a few minutes but happens one time only
try % Load cell parameters to "cellData" structure
  cellData = loadCellParams(cellParam);
catch
  error('Cell parameter file does not exist!');
end
  
if ~exist(ROMfile,'file') % then, need to create the ROM
  try % Load xRA tuning parameters to "xraData" structure
    xraData = loadXRA(configXRA);
  catch
    error('xRA configuration file does not exist!');
  end
  
  ROM = genROM(cellData,xraData,'HRA');
  save(ROMfile,'ROM');
else % load the already-created ROM
  load(ROMfile,'ROM');
end

%% CCCV sim
simData.Ts = 1;
simData.time = 0:1:3000;
simData.SOC0 = 50;
simData.Iapp = -30*ones(3001,1);
simData.T = 25*ones(3001,1);
simData.Papp = 0*simData.Iapp;
simData.warnOff = 1;
ROMout = simROM(ROM,simData,'outBlend',3.0,4.15,0);
save(ROM_CCCV,'ROMout');

%% CPCV sim
simData.Papp = -115*ones(3001,1);
simData.Iapp = 0*simData.Papp;
ROMout = simROM(ROM,simData,'outBlend',3.0,4.15,1);
save(ROM_CPCV,'ROMout');

%% Plotting
close all
load(ROM_CCCV);

figure(1); plot(ROMout.time,100*ROMout.cellSOC); hold on
figure(2); plot(ROMout.time,ROMout.Vcell); hold on
figure(3); plot(ROMout.time,ROMout.Iapp); hold on
figure(4); plot(ROMout.time,ROMout.Vcell.*ROMout.Iapp); hold on

load(ROM_CPCV);
figure(1); plot(ROMout.time,100*ROMout.cellSOC,'--'); 
  xlabel('Time (s)'); ylabel('SOC (%)'); title('State of charge versus time');
  legend('CC/CV','CP/CV','location','northwest'); grid on
  
figure(2); plot(ROMout.time,ROMout.Vcell,'--'); 
  xlabel('Time (s)'); ylabel('Voltage (V)'); title('Terminal voltage versus time');
  legend('CC/CV','CP/CV','location','northwest'); grid on
  
figure(3); plot(ROMout.time,ROMout.Iapp,'--'); 
  xlabel('Time (s)'); ylabel('Current (A)'); title('Cell current versus time');
  legend('CC/CV','CP/CV','location','northwest'); grid on

figure(4); plot(ROMout.time,ROMout.Vcell.*ROMout.Iapp,'--'); 
  xlabel('Time (s)'); ylabel('Power (W)'); title('Cell power versus time');
  legend('CC/CV','CP/CV','location','northwest'); grid on