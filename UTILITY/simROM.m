% function ROMout = simROM(ROM,simData,method,vmin,vmax,pmode)
% 
% Inputs:
%   ROM     = reduced-order model created with an XRA
%   simData = simulation profile loaded using "loadInput"
%   method  = switch between output/model/non blending ("outBlend",
%             "mdlBlend","nonBlend")
%   Vmin    = optional: minimum permitted cell voltage
%   Vmax    = optional: maximum permitted cell voltage
%   Pmode   = optional: if set to 1, input is power and not current
%
% Output:
%   ROMout  = data structure with results from simulation
%
% This utility function executes a reduced-order model using the
% user-prefered blending method
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

function ROMout = simROM(ROM,simData,method,varargin)
  Vmax = Inf; Vmin = -Inf; Pmode = 0;
  if nargin >= 4, Vmin  = varargin{1}; Vmax = varargin{2}; end
  if nargin >= 5, Pmode = varargin{3}; end

  % Update progress
  fprintf('----------------------------------------------------------\n');
  fprintf('Start to simulate ROM (%s)\n',datestr(now)); %#ok<TNOW1,DATST>

  % Simulate ROM using user defined blending methods
  if strcmpi(method,'outBlend')    
    ROMout = outBlend(simData,ROM,Vmin,Vmax,Pmode);
  elseif strcmpi(method,'mdlBlend')
    ROMout = mdlBlend(simData,ROM,Vmin,Vmax,Pmode);
  elseif strcmpi(method,'nonBlend')
    ROMout = nonBlend(simData,ROM,Vmin,Vmax,Pmode);
  end

  fprintf('Finished: %s\n',datestr(now)); %#ok<TNOW1,DATST>
  fprintf('----------------------------------------------------------\n');
end