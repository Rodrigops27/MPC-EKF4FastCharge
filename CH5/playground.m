% This example implements the EKF test cases of chapter 5. It runs with
% either (1) the charge-neutral UDDS profile, (2) the C/5 constant-current
% discharge, and (3) the charge-depleting dynamic profile comprising
% multiple UDDS profiles. Both model blending and output blending are
% demonstrated.

clearvars; close all; clc;
if(~isdeployed),cd(fileparts(which(mfilename)));end
addpath(genpath('../'));
delete(findall(0,'type','figure','tag','TMWWaitbar')); % close all waitbars

TC = 25; % Temperature of simulations
blend = 'MdlB';

SigmaX0 = diag([ones(1,5) 2e6]); % uncertainty of initial state
SigmaW = 1e2;   % uncertainty of current sensor, state equation
SigmaV = 1e-3;  % uncertainty of voltage sensor, output equation; was 1e-3
SOC0 = 25;
% ROMfile    = 'ROMoutB_CCCV.mat'; % generated ROM
ROMfile    = 'ROM_NMC30_HRA.mat'; % generated ROM
FOMoutFile = 'OUT_NMC30_UDDS.mat'; % FOM simulated output
load(FOMoutFile,'FOMout');

% ROMfile    = 'ROM_NMC30_HRA.mat'; % generated ROM
load(ROMfile,'ROM');
% load(ROMfile,'ROMout');
% ROM = ROMout; clear ROMout;

ekfData = initKF(0.95*SOC0,TC,SigmaX0,SigmaV,SigmaW,blend,ROM);
ekfData.maxWarn = 1000; % keep on simulating, even if things look bad
ekfData.trueSOC0 = SOC0/100;

% reserve storage for EKF results... including voltage and SOC (the "+2").
zkEst = NaN(ekfData.nz+2,length(FOMout.Iapp));
zkBound = zkEst;
gamma = NaN(4,length(FOMout.Iapp));
theT = gamma; theZ = gamma;

% hwait = waitbar(0,'Running EKF...');
for k = 1:length(FOMout.Iapp)
    vk = FOMout.Vcell(k);
    ik = FOMout.Iapp(k);
    Tk = FOMout.T(k); % degC
    [zk,zbk,ekfData,Xind] = iterEKF(vk,ik,Tk,ekfData);
    zkEst(:,k)=zk; zkBound(:,k)=zbk;
    gamma(:,k)=Xind.gamma;
    theT(:,k)=Xind.theT;
    theZ(:,k)=Xind.theZ;

    % if mod(k,100)==0, waitbar(k/(length(FOMout.Iapp)-1),hwait);end
end
% close(hwait);

%---------------------------------------------------
% Uncomment if you wish to save data for post-processing
% fileName = sprintf('EKF_%d_%s.mat',simCase,char(blend));
% save(fileName,'FOMout','ekfData','zkEst','zkBound','gamma','theT','theZ');

%% SOC plots
figure;
trueSOC = ekfData.trueSOC0 - cumtrapz(FOMout.Iapp)*ekfData.Ts/(3600*ekfData.Q);
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

figure;
plot(100*(trueSOC' - zkEst(end,:))); hold on; grid on
set(gca,'colororderindex',1); plot(100*zkBound(end,:),':')
set(gca,'colororderindex',1); plot(-100*zkBound(end,:),':')
title(sprintf('SOC estimation error (percent, %s)',char(blend)));
thesisFormat;

count = sum(abs(trueSOC'-zkEst(end,:))>zkBound(end,:));
if any(isnan(zkEst(end,:)))
    fprintf(' - EKF failed (NaN estimates)\n');
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
% thesisFormat;

figure;
plot(FOMout.Vcell-zkEst(end-1,:)'); hold on; grid on
set(gca,'colororderindex',1); plot(+zkBound(end-1,:)',':')
set(gca,'colororderindex',1); plot(-zkBound(end-1,:)',':')
legend('Error','Bounds');
title(sprintf('Voltage estimation error (%s)',char(blend)));
% thesisFormat;

%% Plots of internal variables
% This makes a *lot* of plots and can crash MATLAB, so it is disabled
% by default.

plotKF(zkEst,zkBound,FOMout,ekfData,sprintf('EKF %s',char(blend)));
