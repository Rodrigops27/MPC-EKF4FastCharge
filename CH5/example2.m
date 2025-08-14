% This example implements the SPKF test cases of chapter 5. It runs with
% either (1) the charge-neutral UDDS profile, (2) the C/5 constant-current
% discharge, and (3) the charge-depleting dynamic profile comprising
% multiple UDDS profiles. Both model blending and output blending are
% demonstrated.
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
delete(findall(0,'type','figure','tag','TMWWaitbar')); % close all waitbars

TC = 25; % Temperature of simulations
for simCase = 1:3
  fprintf('Running simulation case %d\n',simCase);
  for blend = {'MdlB','OutB'}

    % Set up SPKF tuning parameters as well as data files and true initial SOC
    switch simCase
      case 1
        SigmaX0 = diag([ones(1,5) 2e6]); % uncertainty of initial state
        SigmaW = 1e2;   % uncertainty of current sensor, state equation
        SigmaV = 1e-3;  % uncertainty of voltage sensor, output equation 
        ROMfile    = 'ROM_NMC30_HRA.mat'; % generated ROM
        FOMoutFile = 'OUT_NMC30_UDDS.mat'; % FOM simulated output
        SOC0 = 60;
      case 2 
        % SigmaX0 = diag([1e5*ones(1,12) 2e6]); % uncertainty of initial state
        SigmaX0 = diag([ones(1,5) 2e6]); % uncertainty of initial state
        SigmaW = 1e2;   % uncertainty of current sensor, state equation
        % SigmaV = 1e-2;  % uncertainty of voltage sensor, output equation 
        SigmaV = 1e-2;  % uncertainty of voltage sensor, output equation 
        ROMfile    = 'ROM_NMC30_HRA.mat'; % generated ROM
        FOMoutFile = 'OUT_NMC30_CC.mat'; % FOM simulated output
        SOC0 = 100;
      case 3
        SigmaX0 = diag([ones(1,5) 2e6]); % uncertainty of initial state
        SigmaW = 1e2;   % uncertainty of current sensor, state equation
        SigmaV = 1e-2;  % uncertainty of voltage sensor, output equation 
        ROMfile    = 'ROM_NMC30_HRA.mat'; % generated ROM
        FOMoutFile = 'OUT_NMC30_Dynamic.mat'; % FOM simulated output
        SOC0 = 95;
    end

    % Load models
    load(ROMfile,'ROM');
    % Load input current, output voltage, etc.
    load(FOMoutFile,'FOMout');

    spkfData = initKF(0.95*SOC0,TC,SigmaX0,SigmaV,SigmaW,blend,ROM);
    spkfData.maxWarn = 1000; % keep on simulating, even if things look bad
    spkfData.trueSOC0 = SOC0/100;

    % reserve storage for SPKF results... including voltage and SOC (the "+2").
    zkEst = NaN(spkfData.nz+2,length(FOMout.Iapp));
    zkBound = zkEst;
    gamma = NaN(4,length(FOMout.Iapp));
    theT = gamma; theZ = gamma;

    hwait = waitbar(0,'Running SPKF...');
    for k = 1:length(FOMout.Iapp)
      vk = FOMout.Vcell(k);
      ik = FOMout.Iapp(k);
      Tk = FOMout.T(k); % degC
      [zk,zbk,spkfData,Xind] = iterSPKF(vk,ik,Tk,spkfData);
      zkEst(:,k)=zk; zkBound(:,k)=zbk;
      gamma(:,k)=Xind.gamma;
      theT(:,k)=Xind.theT;
      theZ(:,k)=Xind.theZ;

      if mod(k,100)==0, waitbar(k/(length(FOMout.Iapp)-1),hwait); end
    end
    close(hwait);

    % Uncomment if you wish to save data for post-processing
    % fileName = sprintf('SPKF_%d_%s.mat',simCase,char(blend));
    % save(fileName,'FOMout','spkfData','zkEst','zkBound','gamma','theT','theZ');

    %% SOC plots
    figure
    trueSOC = spkfData.trueSOC0 - cumtrapz(FOMout.Iapp)*spkfData.Ts/(3600*spkfData.Q);
    plot(100*trueSOC); hold on; grid on
    SOCmdl = NaN(size(trueSOC));
    for k = 1:length(FOMout.Iapp)
      if any(isnan(theT(:,k))), continue, end
      SOCmdl(k) = ROM.ROMmdls(theT(1,k),theZ(1,k)).SOC * gamma(1,k) + ...
                  ROM.ROMmdls(theT(2,k),theZ(2,k)).SOC * gamma(2,k) + ...
                  ROM.ROMmdls(theT(3,k),theZ(3,k)).SOC * gamma(3,k) + ...
                  ROM.ROMmdls(theT(4,k),theZ(4,k)).SOC * gamma(4,k);
    end
    plot(100*SOCmdl,'k--');
    set(gca,'colororderindex',2); plot(100*zkEst(end,:))
    set(gca,'colororderindex',2); plot(100*(zkEst(end,:)+zkBound(end,:)),':')
    set(gca,'colororderindex',2); plot(100*(zkEst(end,:)-zkBound(end,:)),':')
    legend('True SOC','Blend','SOC estimate','Bounds');
    title(sprintf('SOC estimates (percent, %s)',char(blend)));
    thesisFormat;
    
    figure
    plot(100*(trueSOC' - zkEst(end,:))); hold on; grid on
    set(gca,'colororderindex',1); plot(100*zkBound(end,:),':')
    set(gca,'colororderindex',1); plot(-100*zkBound(end,:),':')
    title(sprintf('SOC estimation error (percent, %s)',char(blend)));
    thesisFormat;
    
    count = sum(abs(trueSOC'-zkEst(end,:))>zkBound(end,:));
    if any(isnan(zkEst(end,:)))
      fprintf(' - SPKF failed (NaN estimates)\n');
    else
      fprintf(' - SOC estimate outside of bounds %g %% of the time (%s)\n',...
        count/length(trueSOC)*100,char(blend));
    end
    
    %% voltage plots
    figure
    plot(FOMout.Vcell); hold on; grid on
    set(gca,'colororderindex',2); plot(zkEst(end-1,:))
    set(gca,'colororderindex',2); plot(zkEst(end-1,:)+zkBound(end-1,:),':')
    set(gca,'colororderindex',2); plot(zkEst(end-1,:)-zkBound(end-1,:),':')
    legend('True voltage','Estimate','Bounds');
    title(sprintf('Voltage and estimates (%s)',char(blend)));
    thesisFormat;
    
    figure
    plot(FOMout.Vcell-zkEst(end-1,:)'); hold on; grid on
    set(gca,'colororderindex',1); plot(+zkBound(end-1,:)',':')
    set(gca,'colororderindex',1); plot(-zkBound(end-1,:)',':')
    legend('Error','Bounds');
    title(sprintf('Voltage estimation error (%s)',char(blend)));
    thesisFormat;

    %% Plots of internal variables
    % This makes a *lot* of plots and can crash MATLAB, so it is disabled
    % by default. 

    % plotKF(zkEst,zkBound,FOMout,spkfData,sprintf('SPKF %s',char(blend)));    
  end
end