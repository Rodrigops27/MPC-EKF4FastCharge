% This code creates a set of ROMs for cells having parameters that are
% similar to the "NMC30" cell. Each set of parameters is somewhat
% randomized so that the cells are not identical. The ROMs are used by
% example6 and example7 to simulate battery packs in PCM and SCM formats.
% During operation, the Nyquist plot of the impedance of the NMC30 cell is
% first plotted, and then the Nyquist plots of the impedances of the
% modified cells are overlaid for comparison. The HRA is executed for each
% modified set of cell parameters to create a ROM for that cell.
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

% Set up the workspace and paths
clearvars; close all; clc;
if(~isdeployed),cd(fileparts(which(mfilename)));end
addpath(genpath('../'));

paramDir = './ROM';
outDir = './ROM/ROM_NMC30_%03d.mat';
numROMs = 9; % The number of ROMs to create

% Load XRA configuration for ROM generation
configXRA = 'defaultHRA.xlsx';  % xRA tuning values
try % Load xRA tuning parameters to "xraData" structure
  xraData = loadXRA(configXRA);  
  xraData.SOC = 0:25:100;
  xraData.Tsamp = 1;
  % uncomment the next line if you wish to see debug plots comparing the
  % ideal and ROM frequency responses
  % xraData.debug = 1; 
catch
  error('xRA configuration file does not exist!');
end

% Load the baseline cell parameters 
cellData = loadCellParams('cellNMC30.xlsx');
SOC = 1;            % Default cell state of charge
T = 273.15 + 25;    % Temperature in K for 25 degC
omega = [0, logspace(-6,12,200)]; s = 1j*omega;  

cellData = evalSetpoint(cellData,s,SOC,T);
cellDataTrial = cellData;

% Evaluate transfer functions to determine cell impedance
[phise_tf,~] = tfPhiseInt(s,[0,3],cellData);
[phie_tf,~]  = tfPhie(s,3,cellData);
Z = -(phise_tf(2,:) - phise_tf(1,:) + phie_tf) + cellData.function.const.Rc();

% Nyquist plot of cell impedance
figure(1); clf; h1 = plot(real(1000*Z),-imag(1000*Z),'k'); hold on
xlabel('Real'); h = ylabel('-Imag'); title('Impedance (milliohms)');
CL = lines(200);

% Make ROMs for slight variations on the sample cell parameters
% Step 1: Randomize some parameter values
% Step 2: Make the ROM files
for theTrial = 0:numROMs-1
  % Start step 1
  fprintf('\nCreating randomized cell %d\n',theTrial);  
  cellDataTrial.common = rmfield(cellDataTrial.common,'s');
  
  RF1 = abs(1+0.01*randn(1));
  cellDataTrial.const.Q = cellData.const.Q*RF1;
  fn_value = char(cellData.function.const.Q);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.const.Q = fn_value;      
  
  RF1 = abs(1+0.01*randn(1));
  cellDataTrial.neg.Rf = cellData.neg.Rf*RF1;
  fn_value = char(cellData.function.neg.Rf);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.neg.Rf = fn_value;      

  RF1 = abs(1+0.1*randn(1));
  cellDataTrial.neg.k0 = cellData.neg.k0*RF1;
  fn_value = char(cellData.function.neg.k0);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.neg.k0 = fn_value;      
  
  RF1 = abs(1+0.1*randn(1));
  cellDataTrial.pos.k0 = cellData.pos.k0*RF1;
  fn_value = char(cellData.function.pos.k0);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.pos.k0 = fn_value;      
  
  RF1 = abs(1+0.01*randn(1));
  cellDataTrial.neg.sigma = cellData.neg.sigma*RF1;
  fn_value = char(cellData.function.neg.sigma);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.neg.sigma = fn_value;      
  
  RF1 = abs(1+0.01*randn(1));
  cellDataTrial.pos.sigma = cellData.pos.sigma*RF1;
  fn_value = char(cellData.function.pos.sigma);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.pos.sigma = fn_value;      
  
  RF1 = abs(1+0.01*randn(1));
  cellDataTrial.neg.kappa = cellData.neg.kappa*RF1;
  fn_value = char(cellData.function.neg.kappa);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.neg.kappa = fn_value;      

  RF1 = abs(1+0.01*randn(1));
  cellDataTrial.sep.kappa = cellData.sep.kappa*RF1;
  fn_value = char(cellData.function.sep.kappa);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.sep.kappa = fn_value;      
  
  RF1 = abs(1+0.01*randn(1));
  cellDataTrial.pos.kappa = cellData.pos.kappa*RF1;
  fn_value = char(cellData.function.pos.kappa);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.pos.kappa = fn_value;      

  RF1=abs(1+0.2*randn(1));
  cellDataTrial.neg.Dsref = cellData.neg.Dsref*RF1;
  cellDataTrial.neg.Ds = cellData.neg.Ds*RF1;
  fn_value = char(cellData.function.neg.Dsref);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.neg.Dsref = fn_value;      
  
  RF1=abs(1+0.2*randn(1));
  cellDataTrial.pos.Ds = cellData.pos.Dsref*RF1;
  cellDataTrial.pos.Ds = cellData.pos.Ds*RF1;
  fn_value = char(cellData.function.pos.Dsref);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.pos.Dsref = fn_value;      
  
  RF1=abs(1+0.03*randn(1));
  cellDataTrial.neg.Cdl = cellData.neg.Cdl*RF1;
  fn_value = char(cellData.function.neg.Cdl);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.neg.Cdl = fn_value;      
    
  RF1=abs(1+0.03*randn(1));
  cellDataTrial.pos.Cdl = cellData.pos.Cdl*(1+0.05*randn(1));
  fn_value = char(cellData.function.pos.Cdl);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.pos.Cdl = fn_value;      

  RF1=abs(1+0.01*randn(1));
  cellDataTrial.neg.qe = cellData.neg.qe*RF1;
  fn_value = char(cellData.function.neg.qe);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.neg.qe = fn_value;      
  
  RF1=abs(1+0.01*randn(1));
  cellDataTrial.sep.qe = cellData.sep.qe*RF1;
  fn_value = char(cellData.function.sep.qe);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.sep.qe = fn_value;      
  
  RF1=abs(1+0.01*randn(1));
  cellDataTrial.pos.qe = cellData.pos.qe*RF1;
  fn_value = char(cellData.function.pos.qe);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.pos.qe = fn_value;      

  % For the Rdl components, I multiply the standard values by two. I find
  % that this reduces the tendency to oscillation in the simulations (which
  % is an artifact of the discrete-time nature of the simulation not being
  % able to capture the implicit differentiation of the double-layer
  % elements--a very high frequency effect--very well). Increasing Rdl
  % dampens this effect.
  RF1=2;
  cellDataTrial.neg.Rdl = cellData.neg.Rdl*RF1;
  fn_value = char(cellData.function.neg.Rdl);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.neg.Rdl = fn_value;      
  
  RF1=2;
  cellDataTrial.pos.Rdl = cellData.pos.Rdl*RF1;
  fn_value = char(cellData.function.pos.Rdl);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      theExpr = fn_value(fn_ind2+1:end);
      fn_value = sprintf('@(x,T)(%g*(%s))',RF1,theExpr);
  cellDataTrial.function.pos.Rdl = fn_value;      

  % Prime the transfer functions with the common computations
  [C,L,J,Z,Rct] = tfCommon(s,cellDataTrial);
  cellDataTrial.common.s = s;
  cellDataTrial.common.C = C;
  cellDataTrial.common.L = L;
  cellDataTrial.common.J = J;
  cellDataTrial.common.Z = Z;
  cellDataTrial.common.Rct = Rct;
  cellDataTrial.common.ind = ['L1n=1;L2n=2;L1s=3;L1p=4;L2p=5;' ...
      'c1n=1;c2n=2;c3n=3;c4n=4;c1s=5;c2s=6;c1p=7;c2p=8;c3p=9;c4p=10;'...
      'j1n=1;j2n=2;j3n=3;j4n=4;j1p=5;j2p=6;j3p=7;j4p=8;' ...
      'Zsen=1;Zsep=2;Zsn=3;Zsp=4;Isn=5;Isp=6;'];

  % Evaluate transfer functions to determine cell impedance
  [phise_tf,~] = tfPhiseInt(s,[0,3],cellDataTrial);
  [phie_tf,~]  = tfPhie(s,3,cellDataTrial);
  ZDOT = -(phise_tf(2,:) - phise_tf(1,:) + phie_tf) + cellDataTrial.function.const.Rc();

  figure(1); 
  plot(1000*real(ZDOT),-imag(1000*ZDOT),'color',CL(theTrial+1,:),'markersize',12); 
  ylim([-0.1 2.1]); grid on; axis square; drawnow

  % fileName = sprintf('ROM/TRIAL%03d.mat',theTrial);
  % save(fileName,'cellDataTrial');

  % Start step 2
  fprintf(' - Generating ROM for randomized cell %d\n',theTrial);
  
  % Update format of "const" from character string to function format
  FF = fields(cellDataTrial.function.const);
  for k = 1:length(FF)
    if ischar(cellDataTrial.function.const.(FF{k}))
      cellDataTrial.function.const.(FF{k}) = eval(cellDataTrial.function.const.(FF{k}));
    end
  end  
  % Update format of "neg" 
  FF = fields(cellDataTrial.function.neg);
  for k = 1:length(FF)
    if ischar(cellDataTrial.function.neg.(FF{k}))
      cellDataTrial.function.neg.(FF{k}) = eval(cellDataTrial.function.neg.(FF{k}));
    end
  end  
  % Update format of "sep" 
  FF = fields(cellDataTrial.function.sep);
  for k = 1:length(FF)
    if ischar(cellDataTrial.function.sep.(FF{k}))
      cellDataTrial.function.sep.(FF{k}) = eval(cellDataTrial.function.sep.(FF{k}));
    end
  end
  % Update format of "pos" 
  FF = fields(cellDataTrial.function.pos);
  for k = 1:length(FF)
    if ischar(cellDataTrial.function.pos.(FF{k}))
      cellDataTrial.function.pos.(FF{k}) = eval(cellDataTrial.function.pos.(FF{k}));
    end
  end
  
  ROM = genROM(cellDataTrial,xraData,'HRA');
  outFile = sprintf(outDir,theTrial);
  save(outFile,'ROM');
end
