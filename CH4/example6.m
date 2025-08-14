% This example reproduces results from section 4.11 of the book. It
% simulates a battery pack built from PCMs. The pack comprises Np cells
% wired in parallel to make PCMs; it then has Ns PCMs wired in series to
% build the overall pack structure.
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

outFile = 'ROM/ROMoutPCM.mat'; % Where to save the simulation results
architecture = 'PCM';          % The pack architecture

Ns = 3; % Number of PCMs wired in series
Np = 3; % Number of cells wired in parallel in every PCM

% The pack simulation is slow since optimizations must be done every time
% step to determine the current split between cells in every PCM. When the
% simulation is performed, we save the results in an output file. If this
% file exists, don't re-run the simulation but instead use the precomputed
% results. If you want to run the simulation, then delete (or rename) the
% output file.
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
if ~exist(outFile,'file')
  % Load sufficient ROMs to make a pack of size Ns by Np
  % If you need more ROMs than are provided in the toolbox, run makeROMs.m
  theCell = 0;
  for ks = 1:Ns
    for kp = 1:Np
      fileName = sprintf('ROM/ROM_NMC30_%03d.mat',theCell);
      if exist(fileName,"file")
        ROM = load(fileName,'ROM');    
        ROMs.ROM(ks,kp) = ROM.ROM; 
        theCell = theCell + 1;
      else
        error('Not enough ROMs exist. Run makeROMs.m to create more.');
      end
    end
  end

  % Load input-current profile used for the pack
  simData = loadInput('inputCyclePack.xlsx');

  % Set up individual initial SOCs for every cell (here, make all the same)
  simData.SOC0 = repmat(simData.SOC0,Ns,Np);

  % Simulate the pack and save the output... this is slow
  ROMout = simROMpack(ROMs,simData,architecture,'outBlend'); 
  save(outFile,'ROMout');
end

%% Plot results
close all
load('ROM/ROMoutPCM.mat','ROMout'); % Load the precomputed simulation data
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle');

[Ns,Np] = size(ROMout); % Battery-pack size
t = ROMout(1,1).time;   % Time vector for plotting

% Plot voltages of all PCMs
for ks=1:Ns
  figure(ks)
  for kp=1:Np
    plot(t/60,ROMout(ks,kp).Vcell); hold on    
  end
  title(sprintf('All cell voltages for PCM %d',ks));
  xlabel('Time (min)'); ylabel('Voltage (V)'); grid on
end

% Plot voltages of all cells
figure(Ns+1); clf
for ks=1:Ns
  for kp=1:Np
    plot(t/60,ROMout(ks,kp).Vcell); hold on
  end
end
title('All cell voltages for all PCMs');
xlabel('Time (min)'); ylabel('Voltage (V)'); grid on

% Plot currents of all PCMs
for ks=1:Ns
  figure(Ns+1+ks); clf
  for kp=1:Np
    plot(t/60,ROMout(ks,kp).Iapp); hold on
  end
  title(sprintf('All cell currents for PCM %d',ks));
  xlabel('Time (min)'); ylabel('Current (A)'); grid on
end

% Plot currents for all cells in all PCMs
figure(2*Ns+2)
for ks=1:Ns
  for kp=1:Np
    plot(t/60,ROMout(ks,kp).Iapp); hold on, 
  end
end
title('All cell currents for all PCMs');
xlabel('Time (min)'); ylabel('Current (A)'); grid on
 
% Side-reaction overpotential drives side reactions such as SEI-layer
% growth and lithium plating. We desire for its value to be above 0V. If it
% falls below 0V, then lithium plating will occur.
% eta_s = phi_s - phi_e - Rf * i_{f+dl}
figure(2*Ns+3); clf
for ks=1:Ns
  for kp=1:Np
    % In these sims, the negPhise and negIfdl variables are recorded at
    % positions 0 and 1; we desire only position 1 (negative-electrode/
    % separator boundary)
    phise = ROMout(ks,kp).negPhise(:,2);
    ifdl  = ROMout(ks,kp).negIfdl(:,2);
    Rf    = ROMout(ks,kp).cellData.neg.Rf;
    etas  = phise - Rf*ifdl;
    
    plot(t/60,1000*etas); hold on
  end
end
title('Side-reaction overpotential for all cells in all PCMs');
xlabel('Time (min)'); ylabel('Overpotential (mV)'); grid on