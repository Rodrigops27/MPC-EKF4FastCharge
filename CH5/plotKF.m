% function rmsTable = plotKF(zkEst,zkBound,FOMout,kfData,name)
% 
% Inputs:
%   zkEst   = estimates of different model variables
%   zkBound = 3-sigma bounds on the model variables
%   FOMout  = the COMSOL simulation "truth" data
%   kfData  = configuration of the xKF 
%   name    = 'MdlB' or 'OutB' (for plot titles) 
%
% This is a utility function that plots internal variable estimates,
% bounds, and estimation errors. Returns RMSE table at 4 locations and 8
% variables.
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

function rmsTable = plotKF(zkEst,zkBound,FOMout,kfData,name)
  
  M = 1; % decimation factor
  if length(zkEst)>10000, M=10; end
  CL = lines; CL(3,:) = []; % line colors; delete yellow (too faint)
    
  t = FOMout.time/60; % [mins]
  locations = {'neg:cc','neg:sep','sep:pos','pos:cc'};
  vars = {'Thetae','Phie','Thetass','Phis','Phise','Ifdl','If'};
  rmsTable = zeros(4,8); 
               
  %% ------------- Find the four locations in FOM data -------------
  FOMind = zeros(4,1);
  KFind = zeros(4,1);

  %% ------------------ Thetae ------------------
  xLocs = FOMout.xLocs.Thetae;
  FOMind(1) = find(xLocs == 0,1); % negative current-collector
  FOMind(2) = find(xLocs == 1,1); % neg/sep boundary
  FOMind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  FOMind(4) = find(xLocs == 3,1); % positive current-collector
  xLocs = kfData.loc.Thetae;
  KFind(1) = kfData.ind.Thetae(find(xLocs == 0,1)); % negative cc
  KFind(2) = kfData.ind.Thetae(find(xLocs == 1,1)); % neg/sep boundary
  KFind(3) = kfData.ind.Thetae(find(xLocs == 2,1)); % pos/sep boundary
  KFind(4) = kfData.ind.Thetae(find(xLocs == 3,1)); % positive cc

  rmsVal = plotFour(FOMout.Thetae(:,FOMind),zkEst(KFind,:)',zkBound(KFind,:)',1,name);
  fprintf('Thetae RMSE: %.10f [ul]\n',rmsVal);
  rmsTable(:,5) = rmsVal;
  
  %% ------------------ Phie ------------------
  xLocs = FOMout.xLocs.Phie;
  FOMind(1) = find(xLocs == 0,1); % negative current-collector
  FOMind(2) = find(xLocs == 1,1); % neg/sep boundary
  FOMind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  FOMind(4) = find(xLocs == 3,1); % positive current-collector
  xLocs = kfData.loc.Phie;
  % 1 = negative phise(0)... handle separately
  KFind(2) = kfData.ind.Phie(find(xLocs == 1,1)); % neg/sep boundary
  KFind(3) = kfData.ind.Phie(find(xLocs == 2,1)); % pos/sep boundary
  KFind(4) = kfData.ind.Phie(find(xLocs == 3,1)); % positive cc

  zkInd = kfData.ind.Phise0;
  ROMdata = [-zkEst(zkInd,:)', zkEst(KFind(2:end),:)'];
  ROMbnds = [zkBound(zkInd,:)', zkBound(KFind(2:end),:)'];
  rmsVal = plotFour(FOMout.Phie(:,FOMind),ROMdata,ROMbnds,2,name);
  fprintf('Phie RMSE: %.10f mV\n',rmsVal*1000);
  rmsTable(:,7) = rmsVal*1000;

  %% ------------------ Thetass ------------------
  xLocs = FOMout.xLocs.Thetass;
  FOMind(1) = find(xLocs == 0,1); % negative current-collector
  FOMind(2) = find(xLocs == 1,1); % neg/sep boundary
  FOMind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  FOMind(4) = find(xLocs == 3,1); % positive current-collector
  xLocs = kfData.loc.Thetass;
  KFind(1) = kfData.ind.Thetass(find(xLocs == 0,1)); % negative current collector (N/A)
  KFind(2) = kfData.ind.Thetass(find(xLocs == 1,1)); % neg/sep boundary
  KFind(3) = kfData.ind.Thetass(find(xLocs == 2,1)); % pos/sep boundary
  KFind(4) = kfData.ind.Thetass(find(xLocs == 3,1)); % positive cc

  rmsVal = plotFour(FOMout.Thetass(:,FOMind),zkEst(KFind,:)',zkBound(KFind,:)',3,name);
  fprintf('Thetass RMSE: %.10f [ul]\n',rmsVal);
  rmsTable(:,4) = rmsVal;
  
  %% ------------------ Phis ------------------
  xLocs = FOMout.xLocs.Phis;
  FOMind(1) = find(xLocs == 0,1); % negative current-collector
  FOMind(2) = find(xLocs == 1,1); % neg/sep boundary
  FOMind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  FOMind(4) = find(xLocs == 3,1); % positive current-collector
  xLocs = kfData.loc.Phis;
  % KFind(1) = find(xLocs == 0,1); % negative current collector (N/A)
  KFind(2) = kfData.ind.Phis(find(xLocs == 1,1)); % neg/sep boundary
  KFind(3) = kfData.ind.Phis(find(xLocs == 2,1)); % pos/sep boundary
  KFind(4) = size(zkEst,1)-1;

  ROMdata = [zeros(size(FOMout.Phis(:,1))), zkEst(KFind(2:end),:)'];
  ROMbnds = [zeros(size(FOMout.Phis(:,1))), zkBound(KFind(2:end),:)'];
  rmsVal = plotFour(FOMout.Phis(:,FOMind),ROMdata,ROMbnds,4,name);
  fprintf('Phis RMSE: %.10f mV\n',rmsVal*1000);
  rmsTable(:,6) = rmsVal*1000;

  %% ------------------ Phise ------------------
  xLocs = FOMout.xLocs.Phise;
  FOMind(1) = find(xLocs == 0,1); % negative current-collector
  FOMind(2) = find(xLocs == 1,1); % neg/sep boundary
  FOMind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  FOMind(4) = find(xLocs == 3,1); % positive current-collector
  xLocs = kfData.loc.Phise;
  KFind(1) = kfData.ind.Phise(find(xLocs == 0,1)); % negative current collector (N/A)
  KFind(2) = kfData.ind.Phise(find(xLocs == 1,1)); % neg/sep boundary
  KFind(3) = kfData.ind.Phise(find(xLocs == 2,1)); % pos/sep boundary
  KFind(4) = kfData.ind.Phise(find(xLocs == 3,1)); % positive cc

  rmsVal = plotFour(FOMout.Phise(:,FOMind),zkEst(KFind,:)',zkBound(KFind,:)',5,name);
  fprintf('Phise RMSE: %.10f mV\n',rmsVal*1000);
  rmsTable(:,8) = rmsVal*1000;

  %% ------------------ Ifdl ------------------
  xLocs = FOMout.xLocs.Ifdl;
  FOMind(1) = find(xLocs == 0,1); % negative current-collector
  FOMind(2) = find(xLocs == 1,1); % neg/sep boundary
  FOMind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  FOMind(4) = find(xLocs == 3,1); % positive current-collector
  xLocs = kfData.loc.Ifdl;
  KFind(1) = kfData.ind.Ifdl(find(xLocs == 0,1)); % negative current collector (N/A)
  KFind(2) = kfData.ind.Ifdl(find(xLocs == 1,1)); % neg/sep boundary
  KFind(3) = kfData.ind.Ifdl(find(xLocs == 2,1)); % pos/sep boundary
  KFind(4) = kfData.ind.Ifdl(find(xLocs == 3,1)); % positive cc

  rmsVal = plotFour(FOMout.Ifdl(:,FOMind),zkEst(KFind,:)',zkBound(KFind,:)',6,name);
  fprintf('Ifdl RMSE: %.10f A\n',rmsVal);
  rmsTable(:,1) = rmsVal;

  %% ------------------ If ------------------
  xLocs = FOMout.xLocs.If;
  FOMind(1) = find(xLocs == 0,1); % negative current-collector
  FOMind(2) = find(xLocs == 1,1); % neg/sep boundary
  FOMind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  FOMind(4) = find(xLocs == 3,1); % positive current-collector
  xLocs = kfData.loc.If;
  KFind(1) = kfData.ind.If(find(xLocs == 0,1)); % negative current collector (N/A)
  KFind(2) = kfData.ind.If(find(xLocs == 1,1)); % neg/sep boundary
  KFind(3) = kfData.ind.If(find(xLocs == 2,1)); % pos/sep boundary
  KFind(4) = kfData.ind.If(find(xLocs == 3,1)); % positive cc

  rmsVal = plotFour(FOMout.If(:,FOMind),zkEst(KFind,:)',zkBound(KFind,:)',7,name);
  fprintf('If RMSE: %.10f A\n',rmsVal);
  rmsTable(:,2) = rmsVal;
  
  %% ------------------ PLOT!! ------------------
  function rmsVal = plotFour(FOMsigs,xKFsigs,xKFbnds,theVar,names)
    rmsVal = zeros(4,1);
    for xx = 1:4
      figure;
      clf; plot(t(1:M:end),FOMsigs(1:M:end,xx),'k-','DisplayName','FOM'); hold on;
      plot(t(1:M:end),xKFsigs(1:M:end,xx),'DisplayName',names,'color',CL(1,:));
      
      error = FOMsigs(:,xx) - xKFsigs(:,xx);
      legend('NumColumns',2);
      xlabel('Time (min)'); ylabel('Value');
      title(sprintf('Estimate for %s at %s',vars{theVar},locations{xx}));
      grid on; xlim([0 t(end)]);
      thesisFormat;
      % Compute RMSE
      rmsVal(xx)   = rms(error);
    end

    for xx = 1:4 % Plot estimation errors
      figure;
      error = FOMsigs(:,xx) - xKFsigs(:,xx);
      set(gca,'colororderindex',1);  
      plot(t(1:M:end),error(1:M:end),'DisplayName',names,'color',CL(1,:)); hold on;
      h1 = plot(t(1:M:end),+xKFbnds(1:M:end,xx),'color',CL(1,:));
      h1.Annotation.LegendInformation.IconDisplayStyle = 'off';      
      h2 = plot(t(1:M:end),-xKFbnds(1:M:end,xx),'color',CL(1,:));
      h2.Annotation.LegendInformation.IconDisplayStyle = 'off';      
      
      legend; xlabel('Time (min)'); ylabel('Value');
      title(sprintf('Error for %s at %s',vars{theVar},locations{xx}));
      grid on; xlim([0 t(end)]);
      % thesisFormat;
      set(h1,'linewidth',1);
      set(h2,'linewidth',1);
    end
  end
end