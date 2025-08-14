% function plotResults(ROMout,FOMout,varargin)
% 
% Inputs:
%   ROMout  = output structure from ROM simulation
%   FOMout  = output structure from FOM simulation
%   simCase = (optional) a flag used for positioning legends
%   M       = (optional) a decimation factor to reduce number of points
%             plotted
%
% This is a utility function that accompanies the "example2" through
% "example4" scripts. It is used to plot the results of FOM and ROM
% simulations for visual comparison. It also displays RMSE between the FOM
% and ROM predictions.
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

function plotResults(ROMout,FOMout,varargin)
  if nargin < 3
    simCase = 1;
  else
    simCase = varargin{1};
  end
  if nargin < 4
    M = 1;
  else
    M = varargin{2};
  end
  colors = {[0,0.4470,0.7410]; % blueish
            [0.4660,0.6740,0.1880]; % greenish
            [0.4940,0.1840,0.5560]}; % purplish

  HWIDTH = 1; % linewidth of FOM result
  MSIZE = 20; % ROM marker size

  tFOM = FOMout.time/60;
  tROM = ROMout.time/60; % [mins]
  locations = {' at neg/cc boundary',' at neg/sep boundary',...
               ' at pos/sep boundary',' at pos/cc boundary'};

  rmsTable = zeros(4,8);

  %% ---------- Input current and output voltage ------------------
  figure;clf;
  plot(tFOM(1:M:end),FOMout.Iapp(1:M:end),'color',colors{1});hold on;

  xlabel('Time (min)');ylabel('Current (A)');
  title('Cell current profile'); grid on;

  % Compute RMSE
  error = FOMout.Iapp(:) - ROMout.Iapp(:);
  rmse  = rms(error);
  fprintf('\nIapp RMSE:    %.10f A\n',rmse);

  figure;clf;
  plot(tFOM(1:M:end),FOMout.Vcell(1:M:end),'color',colors{1});hold on;
  plot(tROM(1:M:end),ROMout.Vcell(1:M:end),'.','markersize',MSIZE,'color',colors{2});
  h = plot(tROM(1:M:end),FOMout.Vcell(1:M:end),'color',colors{1});
  switch simCase
    case 1, legendLoc='northeast';
    case 2, legendLoc='northeast';
    case 3, legendLoc='southeast';
  end
  legend('FOM','ROM','location',legendLoc);

  xlabel('Time (min)');ylabel('Voltage (V)');
  title(sprintf('Cell terminal voltage (%s)',ROMout.blending));
  grid on; set(h,'linewidth',HWIDTH);

  % Compute RMSE
  error = FOMout.Vcell(:) - ROMout.Vcell(:);
  vrms   = rms(error);
  fprintf('Vcell RMSE:   %.10f mV\n',vrms*1e3);

  
  %% ------------- Find the four locations in FOM data -------------
  ind = zeros(4,1);

  %% ------------------ Phis ------------------
  xLocs = FOMout.xLocs.Phis;
  ind(1) = find(xLocs == 0,1); % negative current-collector
  ind(2) = find(xLocs == 1,1); % neg/sep boundary
  ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  ind(4) = find(xLocs == 3,1); % positive current-collector

  for x = 1:2
    figure;clf;
    plot(tFOM(1:M:end),FOMout.Phis(1:M:end,ind(x+1)),'color',colors{1}); hold on;
    plot(tROM(1:M:end),ROMout.Phis(1:M:end,x),'.','markersize',MSIZE,'color',colors{2});
    h = plot(tFOM(1:M:end),FOMout.Phis(1:M:end,ind(x+1)),'color',colors{1});    
    switch simCase
      case 1, if x == 1, legendLoc='southeast'; else, legendLoc='northeast'; end
      case 2, if x == 1, legendLoc='southeast'; else, legendLoc='northeast'; end
      case 3, legendLoc='southeast';
    end        
    legend('FOM','ROM','location',legendLoc);
    
    xlabel('Time (min)');ylabel('Potential (V)');
    title(sprintf('Phis %s (%s)',locations{x+1},ROMout.blending));
    grid on; set(h,'linewidth',HWIDTH);  

    error = FOMout.Phis(1:M:end,ind(x+1)) - ROMout.Phis(1:M:end,x);
    rmse  = rms(error);
    fprintf('Phis RMSE:    %.10f mV\n',rmse*1000);
    rmsTable(x+1,6) = rmse*1000;
  end
  rmsTable(4,6) = vrms*1000;

 
  %% ------------------ Phie ------------------
  xLocs = FOMout.xLocs.Phie;
  ind(1) = find(xLocs == 0,1); % negative current-collector
  ind(2) = find(xLocs == 1,1); % neg/sep boundary
  ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  ind(4) = find(xLocs == 3,1); % positive current-collector

  for x = 1:4
    figure;clf;
    plot(tFOM(1:M:end),FOMout.Phie(1:M:end,ind(x)),'color',colors{1}); hold on;
    plot(tROM(1:M:end),ROMout.Phie(1:M:end,x),'.','markersize',MSIZE,'color',colors{2});
    h = plot(tROM(1:M:end),FOMout.Phie(1:M:end,ind(x)),'color',colors{1});        
    switch simCase
      case 1, legendLoc='northeast';
      case 2, legendLoc='northeast';
      case 3, legendLoc='southeast';
    end                
    legend('FOM','ROM','location',legendLoc);
    
    xlabel('Time (min)');ylabel('Potential (V)');
    title(sprintf('Phie %s (%s)',locations{x},ROMout.blending));
    grid on; set(h,'linewidth',HWIDTH);  

    error = FOMout.Phie(1:M:end,ind(x)) - ROMout.Phie(1:M:end,x);
    rmse  = rms(error);
    fprintf('Phie RMSE:    %.10f mV\n',rmse*1000);
    rmsTable(x,7) = rmse*1000;
  end

  %% ------------------ Phise ------------------
  xLocs = FOMout.xLocs.Phise;
  ind(1) = find(xLocs == 0,1); % negative current-collector
  ind(2) = find(xLocs == 1,1); % neg/sep boundary
  ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  ind(4) = find(xLocs == 3,1); % positive current-collector

  for x = 1:4
    figure;clf;
    plot(tFOM(1:M:end),FOMout.Phise(1:M:end,ind(x)),'color',colors{1}); hold on;
    plot(tROM(1:M:end),ROMout.Phise(1:M:end,x),'.','markersize',MSIZE,'color',colors{2});
    h = plot(tFOM(1:M:end),FOMout.Phise(1:M:end,ind(x)),'color',colors{1});        
    switch simCase
      case 1, if x>2, legendLoc='northeast';
        else, legendLoc='southeast'; end
      case 2, if x>2, legendLoc='northeast';
        else, legendLoc='southeast'; end
      case 3, if x<3, legendLoc='northeast';
        else, legendLoc='southeast'; end
    end        
    legend('FOM','ROM','location',legendLoc);
    
    xlabel('Time (min)');ylabel('Potential (V)');
    title(sprintf('Phise %s (%s)',locations{x},ROMout.blending));
    grid on; set(h,'linewidth',HWIDTH);  

    error = FOMout.Phise(1:M:end,ind(x)) - ROMout.Phise(1:M:end,x);
    rmse  = rms(error);
    fprintf('Phise RMSE:   %.10f mV\n',rmse*1000);
    rmsTable(x,8) = rmse*1000;
  end

  %% ------------------ Thetae ------------------
  xLocs = FOMout.xLocs.Thetae;
  ind(1) = find(xLocs == 0,1); % negative current-collector
  ind(2) = find(xLocs == 1,1); % neg/sep boundary
  ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  ind(4) = find(xLocs == 3,1); % positive current-collector

  for x = 1:4
    figure;clf;
    plot(tFOM(1:M:end),FOMout.Thetae(1:M:end,ind(x)),'color',colors{1}); hold on;
    plot(tROM(1:M:end),ROMout.Thetae(1:M:end,x),'.','markersize',MSIZE,'color',colors{2});
    h = plot(tFOM(1:M:end),FOMout.Thetae(1:M:end,ind(x)),'color',colors{1});            
    switch simCase
      case 1, if x>2, legendLoc='northwest';
        else, legendLoc='southeast'; end
      case 2, if x>2, legendLoc='northeast';
        else, legendLoc='southeast'; end
      case 3, if x<3, legendLoc='northeast';
        else, legendLoc='southeast'; end
    end        
    legend('FOM','ROM','location',legendLoc);
    
    xlabel('Time (min)');ylabel('Concentration ratio (u/l)');
    title(sprintf('Thetae %s (%s)',locations{x},ROMout.blending));
    grid on; set(h,'linewidth',HWIDTH);  

    error = FOMout.Thetae(1:M:end,ind(x)) - ROMout.Thetae(1:M:end,x);
    rmse  = rms(error);
    fprintf('Thetae RMSE:  %.10f [unitless]\n',rmse);
    rmsTable(x,5) = rmse;
  end

  %% ------------------ Thetass ------------------
  xLocs = FOMout.xLocs.Thetass;
  ind(1) = find(xLocs == 0,1); % negative current-collector
  ind(2) = find(xLocs == 1,1); % neg/sep boundary
  ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  ind(4) = find(xLocs == 3,1); % positive current-collector

  for x = 1:4
    figure;clf;
    plot(tFOM(1:M:end),FOMout.Thetass(1:M:end,ind(x)),'color',colors{1}); hold on;
    plot(tROM(1:M:end),ROMout.Thetass(1:M:end,x),'.','markersize',MSIZE,'color',colors{2});
    h = plot(tFOM(1:M:end),FOMout.Thetass(1:M:end,ind(x)),'color',colors{1});            
    switch simCase
      case 1, if x>2, legendLoc='southeast';
        else, legendLoc='northeast'; end
      case 2, if x>2, legendLoc='southeast';
        else, legendLoc='northeast'; end
      case 3, if x<3, legendLoc='southeast';
        else, legendLoc='northeast'; end
    end        
    legend('FOM','ROM','location',legendLoc);
    
    xlabel('Time (min)');ylabel('Concentration ratio (u/l)');
    title(sprintf('Thetass %s (%s)',locations{x},ROMout.blending));
    grid on; set(h,'linewidth',HWIDTH);  

    error = FOMout.Thetass(1:M:end,ind(x)) - ROMout.Thetass(1:M:end,x);
    rmse  = rms(error);
    fprintf('Thetass RMSE: %.10f [unitless]\n',rmse);
    rmsTable(x,4) = rmse;
  end

  %% ------------------ Ifdl ------------------
  xLocs = FOMout.xLocs.Ifdl;
  ind(1) = find(xLocs == 0,1); % negative current-collector
  ind(2) = find(xLocs == 1,1); % neg/sep boundary
  ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  ind(4) = find(xLocs == 3,1); % positive current-collector

  for x = 1:4
    figure;clf;
    plot(tFOM(1:M:end),FOMout.Ifdl(1:M:end,ind(x)),'color',colors{1}); hold on;
    plot(tROM(1:M:end),ROMout.Ifdl(1:M:end,x),'.','markersize',MSIZE,'color',colors{2});
    h = plot(tFOM(1:M:end),FOMout.Ifdl(1:M:end,ind(x)),'color',colors{1});                
    switch simCase
      case 1, legendLoc='southeast';
      case 2, if x == 2, legendLoc='northeast';
        else, legendLoc='southeast'; end
      case 3, if x<3, legendLoc='northeast';
        else, legendLoc='southeast'; end
    end        
    legend('FOM','ROM','location',legendLoc);
    
    xlabel('Time (min)');ylabel('Interfacial lithium flux (A)');
    title(sprintf('Ifdl %s (%s)',locations{x},ROMout.blending));
    grid on; set(h,'linewidth',HWIDTH);  

    % Compute RMS
    error = FOMout.Ifdl(1:M:end,ind(x)) - ROMout.Ifdl(1:M:end,x);
    rmse  = rms(error);
    fprintf('Ifdl RMSE:    %.10f A\n',rmse);
    rmsTable(x,1) = rmse;
  end

  %% ------------------ If ------------------
  xLocs = FOMout.xLocs.If;
  ind(1) = find(xLocs == 0,1); % negative current-collector
  ind(2) = find(xLocs == 1,1); % neg/sep boundary
  ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  ind(4) = find(xLocs == 3,1); % positive current-collector

  for x = 1:4
    figure;clf;
    plot(tFOM(1:M:end),FOMout.If(1:M:end,ind(x)),'color',colors{1}); hold on;
    plot(tROM(1:M:end),ROMout.If(1:M:end,x),'.','markersize',MSIZE,'color',colors{2});
    h = plot(tFOM(1:M:end),FOMout.If(1:M:end,ind(x)),'color',colors{1});                    
    switch simCase
      case 1, legendLoc='southeast';
      case 2, if x == 2, legendLoc='northeast';
        else, legendLoc='southeast'; end
      case 3, if x<3, legendLoc='northeast';
        else, legendLoc='southeast'; end
    end        
    legend('FOM','ROM','location',legendLoc);
    
    xlabel('Time (min)');ylabel('Interfacial lithium flux (A)');
    title(sprintf('If %s (%s)',locations{x},ROMout.blending));
    grid on; set(h,'linewidth',HWIDTH);  

    % Compute RMS
    error = FOMout.If(1:M:end,ind(x)) - ROMout.If(1:M:end,x);
    rmse  = rms(error);
    fprintf('If RMSE:      %.10f A\n',rmse);
    rmsTable(x,2) = rmse;
  end

  %% ------------------ Idl ------------------
  xLocs = FOMout.xLocs.Idl;
  ind(1) = find(xLocs == 0,1); % negative current-collector
  ind(2) = find(xLocs == 1,1); % neg/sep boundary
  ind(3) = find(xLocs == 2,1,'last'); % pos/sep boundary
  ind(4) = find(xLocs == 3,1); % positive current-collector

  % ROMout.Idl = ROMout.Ifdl - ROMout.If;
  for x = 1:4
    figure;clf;
    plot(tFOM(1:M:end),FOMout.Idl(1:M:end,ind(x)),'color',colors{1}); hold on;
    plot(tROM(1:M:end),ROMout.Idl(1:M:end,x),'.','markersize',MSIZE,'color',colors{2});
    h = plot(tFOM(1:M:end),FOMout.Idl(1:M:end,ind(x)),'color',colors{1});                        
    switch simCase
      case 1, legendLoc='southeast';
      case 2, if x < 3, legendLoc='northeast';
        else, legendLoc='southeast'; end
      case 3, if x<3, legendLoc='northeast';
        else, legendLoc='southeast'; end
    end        
    legend('FOM','ROM','location',legendLoc);
    
    xlabel('Time (min)');ylabel('Interfacial lithium flux (A)');
    title(sprintf('Idl %s (%s)',locations{x},ROMout.blending));
    grid on; set(h,'linewidth',HWIDTH);  

    % Compute RMS
    error = FOMout.Idl(1:M:end,ind(x)) - ROMout.Idl(1:M:end,x);
    rmse  = rms(error);
    fprintf('Idl RMSE:     %.10f A\n',rmse);
    rmsTable(x,3) = rmse;
  end
end