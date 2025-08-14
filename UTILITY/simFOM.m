% function [FOM,FOMout] = simFOM(FOM,simData)
% 
% Inputs:
%   FOM     = COMSOL object containing the full-order model, created by
%             genFOM.m
%   simData = simulation profile loaded using "loadInput"
%
% Output:
%   FOM     = Updated COMSOL object containing the full-order model and all
%             simulated data  
%   FOMout  = Data structure with results from simulation
%
% This utility function executes a COMSOL model for specified input-current
% and input-temperature profiles, starting at an initial cell SOC, using 
% the LiveLink for MATLAB interface. The returned COMSOL object can be
% saved to a file using mphsave.m, or loaded into the COMSOL GUI using
% mphlaunch.m. 
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

function [FOM,FOMout] = simFOM(FOM,simData)
  import com.comsol.model.*      %#ok<NSTIMP>
  import com.comsol.model.util.* %#ok<NSTIMP>
  debugFlag = 1; % set to "1" if you want to see progress update, else "0"
  FOMout = []; 
  
  % Shift time vector when COMSOL samples all output values by TSHIFT to
  % allow COMSOL to respond to abrupt changes in input current. 
  %
  % If TSHIFT=0, then COMSOL output does not have time to respond to step-
  % like instant changes in input current before sampling voltage (etc.), 
  % so the voltage change does not appear until the following time sample,
  % making the overall voltage profile appear delayed by one time sample.
  % That is, a nonzero TSHIFT is needed to capture feedthrough resistance
  % (ESR) effects. A TSHIFT of 0.1% of the sample period works; a shift of
  % 1% might be slightly better.
  if isfield(simData,'TSHIFT')
    TSHIFT = simData.TSHIFT;
  else
%     TSHIFT = 1e-2; 
    TSHIFT = 0.0025;
  end
                 
  updateInputs;
  runStudy;
  collectResults;
  msg('\n');
  
  function updateInputs
    msg('Updating inputs in FOM...');
    tvec = simData.time(:);
    ivec = simData.Iapp(:);
    Tvec = simData.T(:);
    
    % Replace default input-current function
    if length(tvec) == length(ivec),  ivec = ivec(1:end-1); end
    v1 = tvec(1:end-1); v2 = tvec(2:end); v3 = ivec;
    Ides = sprintf('''%g'' ''%g'' ''%g'';',[v1(:),v2(:),v3(:)]');
    Ides = eval(sprintf('{ %s }',Ides));
    FOM.func('Ides').set('pieces', Ides);
    
    % Replace default input temperature function
    if length(tvec) == length(Tvec),  Tvec = Tvec(1:end-1); end
    v1 = tvec(1:end-1); v2 = tvec(2:end); v3 = Tvec+273.15;
    Temperature = sprintf('''%g'' ''%g'' ''%g'';',[v1(:),v2(:),v3(:)]');
    Temperature = eval(sprintf('{ %s }',Temperature));
    FOM.func('Temperature').set('pieces', Temperature);

    % Replace default initial SOC
    FOM.param.set('z0', num2str(simData.SOC0/100), 'Initial cell SOC');
    
    % Replace default study values
    tvar = tvec+TSHIFT; 
    FOM.study('std1').feature('time').set('tlist', tvar);      
    FOM.sol('sol1').feature('t1').set('tlist', tvar);
  end   % updateInputs
  function runStudy 
    import com.comsol.model.* 
    import com.comsol.model.util.*    
    msg('\nRunning study...');
    if ~ismac, ModelUtil.showProgress(true); end    
    % mphlaunch(mphtags(FOM));    
    try
      FOM.sol('sol1').runAll;
    catch
      warning('COMSOL did not converge! Results computed will be saved.')
    end    
  end      % runStudy
  function collectResults
    msg('\nCollecting results...');    
    % In negative electrode and at its current collector
    vars = {'phi_s','phi_e','theta_e','if','idl','thetass','ifdl'};
    data_neg = mpheval(FOM,vars,'Edim',1,'Selection',1);
    locs_neg = data_neg.p(:);
    FOMout.negIfdl    = data_neg.d7;
    FOMout.negIf      = data_neg.d4;
    FOMout.negIdl     = data_neg.d5;
    FOMout.negPhis    = data_neg.d1;
    FOMout.negPhise   = data_neg.d1 - data_neg.d2;
    FOMout.negThetass = data_neg.d6;

    FOMout.negIfdl0    = FOMout.negIfdl(:,1); 
    FOMout.negIf0      = FOMout.negIf(:,1);
    FOMout.negIdl0     = FOMout.negIdl(:,1);
    FOMout.negPhise0   = FOMout.negPhise(:,1);
    FOMout.negThetass0 = FOMout.negThetass(:,1);
        
    % In the electrolyte
    data_elect = mpheval(FOM,{'theta_e','phi_e'},'Edim',1);
    locs_elect = data_elect.p(:);
    FOMout.Phie   = data_elect.d2;
    FOMout.Thetae = data_elect.d1;
    FOMout.xLocs.Phie   = locs_elect;
    FOMout.xLocs.Thetae = locs_elect;
    
    % In positive electrode and at its current collector
    data_pos = mpheval(FOM,vars,'Edim',1,'Selection',3);
    locs_pos = data_pos.p(:);
    FOMout.posIfdl    = data_pos.d7;
    FOMout.posIf      = data_pos.d4;
    FOMout.posIdl     = data_pos.d5;
    FOMout.posPhis    = data_pos.d1;
    FOMout.posPhise   = data_pos.d1 - data_pos.d2;
    FOMout.posThetass = data_pos.d6;

    FOMout.posIfdl3    = FOMout.posIfdl(:,end); 
    FOMout.posIf3      = FOMout.posIf(:,end); % GLP: 01/04/23
    FOMout.posIdl3     = FOMout.posIdl(:,end);
    FOMout.posPhise3   = FOMout.posPhise(:,end);
    FOMout.posThetass3 = FOMout.posThetass(:,end);

    FOMout.xLocs.Ifdl    = [locs_neg; locs_pos];
    FOMout.xLocs.If      = [locs_neg; locs_pos];
    FOMout.xLocs.Idl     = [locs_neg; locs_pos];
    FOMout.xLocs.Phis    = [locs_neg; locs_pos];
    FOMout.xLocs.Phise   = [locs_neg; locs_pos];
    FOMout.xLocs.Thetass = [locs_neg; locs_pos];
                          
    % Simulation control and primary output
    temp = mpheval(FOM,'T');
    FOMout.T = temp.d1(:,end)-273.15; % Degrees Celsius
    tvar = mpheval(FOM,'t');  tvar = tvar.d1(:,1);
    FOMout.time = tvar-TSHIFT;
    input = mpheval(FOM,'Ides');   
    FOMout.Ides = input.d1(:,end);
    input = mpheval(FOM,'Iapp');   
    FOMout.Iapp = input.d1(:,end);
    input = mpheval(FOM,'Vcell','Edim',0);
    FOMout.Vcell = input.d1(:,end);
    
    % Other cell and electrode quantities
    % ROMout.cellSOC = zeros(duration,1); 
    data_neg = mpheval(FOM,'thetasavg_neg');
    FOMout.negSOC = data_neg.d1;
    data_pos = mpheval(FOM,'thetasavg_pos');
    FOMout.posSOC = data_pos.d1;
    
    % Save each variable
    FOMout.Ifdl    = [FOMout.negIfdl FOMout.posIfdl];
    FOMout.If      = [FOMout.negIf FOMout.posIf];
    FOMout.Idl     = [FOMout.negIdl FOMout.posIdl];
    FOMout.Phis    = [FOMout.negPhis FOMout.posPhis];
    FOMout.Phise   = [FOMout.negPhise FOMout.posPhise];
    FOMout.Thetass = [FOMout.negThetass FOMout.posThetass];
  end % collectResults
  function msg(theText)
    if debugFlag
      fprintf(theText);
    end
  end % msg
end