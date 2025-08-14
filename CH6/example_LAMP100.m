% This example reproduces one of the results in section 6.4 of the book.
% Specifically, it simulates the loss of active material in the positive
% electrode, at 100% SOC.
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

scenario = 'LAMP100';

%% Initialize
Vmax = 4.2;
Vmin = 3.0;

load OCP.mat
theta0n = 0.0299; 
theta100n = 0.8926;
theta0p = 0.8510;
theta100p = 0.1361;
QAh = 10;

Qn = QAh/abs(theta100n - theta0n); % Electrode capacities; constant for LLI
Qp = QAh/abs(theta100p - theta0p);

theta = [theta0n; theta100n; theta0p; theta100p];
Qcell = QAh;
QpEps = Qp;

%% Loop over cycles to calculate values of aging variables
qEps = 0.002; % fractional material loss per iteration
for k = 1:300 % Loop over time periods
  Qp = Qp - qEps*Qp;
  
  deltaQ = 0:1e-3:QAh+0.1; % search over this vector of dis/charge depth

  % "Discharge" to find stoichiometry for the V closest to Vmin
  thetaNeg = theta100n - deltaQ/Qn;
  ocpNeg = interp1(OCPneg.z,OCPneg.v,thetaNeg,'linear','extrap');
  thetaPos = theta100p + deltaQ/Qp;
  ocpPos = interp1(OCPpos.z,OCPpos.v,thetaPos,'linear','extrap');
  ocvCell = ocpPos - ocpNeg;
  theta0n = interp1(ocvCell,thetaNeg,Vmin,'linear','extrap');
  theta0p = interp1(ocvCell,thetaPos,Vmin,'linear','extrap');

  % "Charge" to find stoichiometry for the V closest to Vmax
  thetaNeg = theta0n + deltaQ/Qn;
  ocpNeg = interp1(OCPneg.z,OCPneg.v,thetaNeg,'linear','extrap');
  thetaPos = theta0p - deltaQ/Qp;
  ocpPos = interp1(OCPpos.z,OCPpos.v,thetaPos,'linear','extrap');
  ocvCell = ocpPos - ocpNeg;
  theta100n = interp1(ocvCell,thetaNeg,Vmax,'linear','extrap');
  theta100p = interp1(ocvCell,thetaPos,Vmax,'linear','extrap');
  
  theta = [theta, [theta0n; theta100n; theta0p; theta100p]]; %#ok<*AGROW>
  Qcell = [Qcell, Qn*abs(theta100n - theta0n)];
  QAh = Qcell(end);
  QpEps = [QpEps; Qp];
end

%% Prepare to plot results
% Trim to an 80% maximum total-capacity loss
ind = find(Qcell < 0.8*Qcell(1),1,'first');
theta = theta(:,1:ind);
Qcell = Qcell(1:ind);
QpEps = QpEps(1:ind);
CL = lines; MS = 10;
loss = (Qcell(1) - Qcell)/Qcell(1)*100;

%% Plot marker dots on the negative
figure(2); clf; plot(OCPneg.z,OCPneg.v,'k'); hold on;
vHist000 = interp1(OCPneg.z,OCPneg.v,theta(1,:));
plot(theta(1,:),vHist000,'.','markersize',MS,'color',CL(1,:));
vHist100 = interp1(OCPneg.z,OCPneg.v,theta(2,:));
plot(theta(2,:),vHist100,'.','markersize',MS,'color',CL(2,:));
ylim([0 3]); xlim([0 1]); grid on; 
legend('OCP','theta0n','theta100n','numcolumns',3);
title(sprintf('Movement of negative-electrode boundaries for %s case',scenario));
xlabel('Negative-electrode stoichiometry'); ylabel('Open-circuit potential (V)');

% Arrow for theta0n
offset = -0.06; len = 0.1;
X0 = theta(1,round(ind/2)); Y0 = vHist000(round(ind/2)); [X0f,Y0f] = ds2nfu(X0,Y0);
XM = theta(1,round(0.6*ind)); YM = vHist000(round(0.6*ind)); [XMf,YMf] = ds2nfu(XM,YM);
XN = theta(1,end); YN = vHist000(end);                   [XNf,YNf] = ds2nfu(XN,YN);
psi = atan2(YNf-Y0f,XNf-X0f); % find angle of line
XMf = XMf + offset*cos(psi-pi/2); YMf = YMf + offset*sin(psi-pi/2); % center
X0 = XMf - len*cos(psi)/2; XN = XMf + len*cos(psi)/2;
Y0 = YMf - len*sin(psi)/2; YN = YMf + len*sin(psi)/2;
annotation('arrow',[X0 XN],[Y0 YN],'linewidth',2);

%% Plot marker dots on the positive
figure(3); clf; plot(OCPpos.z,OCPpos.v,'k'); hold on;
vHist000 = interp1(OCPpos.z,OCPpos.v,theta(3,:));
plot(theta(3,:),vHist000,'.','markersize',MS,'color',CL(3,:));
vHist100 = interp1(OCPpos.z,OCPpos.v,theta(4,:));
plot(theta(4,:),vHist100,'.','markersize',MS,'color',CL(4,:));
ylim([3.4 4.4]); xlim([0 1]); grid on;
legend('OCP','theta0p','theta100p','numcolumns',3);
title(sprintf('Movement of positive-electrode boundaries for %s case',scenario));
xlabel('Positive-electrode stoichiometry'); ylabel('Open-circuit potential (V)');

% Arrow for theta0p
offset = -0.03; len = 0.1;
X0 = theta(3,1); Y0 = vHist000(1);                       [X0f,Y0f] = ds2nfu(X0,Y0);
XM = theta(3,round(0.25*ind)); YM = vHist000(round(0.25*ind)); [XMf,YMf] = ds2nfu(XM,YM);
XN = theta(3,round(0.4*ind)); YN = vHist000(round(0.4*ind));                   [XNf,YNf] = ds2nfu(XN,YN);
psi = atan2(YNf-Y0f,XNf-X0f); % find angle of line
XMf = XMf + offset*cos(psi-pi/2); YMf = YMf + offset*sin(psi-pi/2); % center
X0 = XMf - len*cos(psi)/2; XN = XMf + len*cos(psi)/2;
Y0 = YMf - len*sin(psi)/2; YN = YMf + len*sin(psi)/2;
annotation('arrow',[X0 XN],[Y0 YN],'linewidth',2);

%% Plot OCV at different levels of capacity loss
figure(4); clf; z = 0:0.001:1; 
snapshots = [0 5 10 15 20];
for k = 1:length(snapshots)
  ind = find(loss > snapshots(k),1,'first'); % for 100% SOC
  thetas = theta(:,ind);
  thetaNeg = thetas(1) + z*(thetas(2)-thetas(1));
  thetaPos = thetas(3) + z*(thetas(4)-thetas(3));
  ocpNeg = interp1(OCPneg.z,OCPneg.v,thetaNeg,'linear','extrap');
  ocpPos = interp1(OCPpos.z,OCPpos.v,thetaPos,'linear','extrap');
  ocvCell = ocpPos - ocpNeg;
  plot(100*z,ocvCell); hold on
end
legendCell = cellstr(num2str(snapshots', 'Loss = %-d%%'));
legend(legendCell,'location','southeast');
ylim([Vmin Vmax]); 
ylabel('Cell OCV (V)');
xlabel('Cell SOC (%)');
title(sprintf('Cell OCV for %s case',scenario));
grid on

% Plot dOCV at different levels of capacity loss
figure(5); clf;
for k = 1:length(snapshots)
  ind = find(loss > snapshots(k),1,'first'); % for 100% SOC
  thetas = theta(:,ind);
  thetaNeg = thetas(1) + z*(thetas(2)-thetas(1));
  thetaPos = thetas(3) + z*(thetas(4)-thetas(3));
  docpNeg = interp1(OCPneg.z,OCPneg.dv,thetaNeg,'linear','extrap');
  docpPos = interp1(OCPpos.z,OCPpos.dv,thetaPos,'linear','extrap');
  docvCell = docpPos*(thetas(4) - thetas(3)) - docpNeg*(thetas(2) - thetas(1));
  plot(100*z,docvCell); hold on
end
legendCell = cellstr(num2str(snapshots', 'Loss = %-d%%'));
legend(legendCell,'location','northeast');
ylim([0 10]);
ylabel('Cell dOCV/dz (V)');
xlabel('Cell SOC (%)');
title(sprintf('Cell differential OCV for %s case',scenario));
grid on
