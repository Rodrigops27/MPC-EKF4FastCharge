% function ROMout = simROMpack(ROMs,simData,architecture,method)
% 
% Inputs:
%   ROM          = matrix of reduced-order models created with an XRA
%   simData      = simulation profile loaded using "loadInput"
%   architecture = either 'SCM' or 'PCM'; the pack configuraion
%   method       = in the future, switch between output/model/non blending
%                  ("outBlend", "mdlBlend","nonBlend"). Presently, only
%                  "outBlend" is implemented.
%
% Output:
%   ROMout       = matrix of data structures with results from simulation
%
% This utility function executes a Ns by Np battery pack using
% reduced-order model and the user-prefered blending method
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

function ROMout = simROMpack(ROMs,simData,architecture,method)

  % Update progress
  fprintf('----------------------------------------------------------\n');
  fprintf('Start to simulate ROM-based pack using %s (%s)\n',...
    architecture,datestr(now)); %#ok<TNOW1,DATST>

  switch upper(method)
    case 'OUTBLEND'      
      ROMout = outBlendPack(ROMs,simData,architecture);
    otherwise 
      error('simROMpack presently works only with output blending.')
  end

  fprintf('Finished: %s\n',datestr(now)); %#ok<TNOW1,DATST>
  fprintf('----------------------------------------------------------\n');
end