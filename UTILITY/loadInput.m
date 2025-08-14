% function simData = loadInput(fileName)
% 
% Inputs:
%   fileName = filename of Excel spreadsheet to be loaded
%
% Outputs:
%   simData  = profile of current and temperature versus time; init SOC
%
% This function loads an Excel spreadsheet that specifies current and
% temperature versus time (plus initial SOC) for a ROM/FOM simulation.
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

function simData = loadInput(fileName)
  [num,~,~] = xlsread(fileName); %#ok<XLSRD>
  simData.time = num(:,1);
  simData.time(isnan(num(:,1))) = [];
  simData.Iapp = num(:,2);
  simData.Iapp(isnan(num(:,2))) = [];
  simData.T    = num(:,3);
  simData.T(isnan(num(:,3))) = [];
  simData.Ts   = abs(num(2,1)-num(1,1));
  simData.SOC0 = num(:,4:5);
  simData.SOC0(isnan(num(:,4:5))) = [];
  simData.SOC0 = simData.SOC0(1); % use only first number, if multiple
end
