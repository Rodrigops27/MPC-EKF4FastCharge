% function ROM = genROM(cellData,xraData,method)
% 
% Inputs:
%   cellData = cell-parameter data structure loaded with "loadCellParams"
%   xraData  = xra-control parameters loaded with "loadXRA"
%   method   = xra method (e.g., 'HRA')
%
% Outputs:
%   ROM      = a structure containing the ROM data at multiple setpoints
%
% This function creates a reduced-order model for the cell defined by
% "cellData" using the XRA tuning parameters in "xraData" using the
% realization algorithm defined by "method"
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

function ROM = genROM(cellData,xraData,method)
  % Build an output data structure has total (#T times #SOC) cells
  out = cell(length(xraData.T),length(xraData.SOC));

  % Loop over every temperature (tt) and SOC (zz) setpoint
  for tt = 1:length(xraData.T)
    for zz = 1:length(xraData.SOC)
      T   = xraData.T(tt)+273.15; % [°K]
      SOC = xraData.SOC(zz)/100;  % [ul]

      % Execute XRA
      fprintf('Starting: %s\n',datestr(now)); %#ok<TNOW1,DATST>
      fprintf(' - XRA: %s, T: %g°C, SOC: %g%%, Ts: %gs, Order: %g\n',...
          method,T-273.15,SOC*100,xraData.Tsamp,xraData.n);

      methodFile = sprintf('xra%s.m',upper(method));
      if ~exist(methodFile,'file')
        error('Invalid XRA method given');
      end
      execXRA = sprintf('[out{tt,zz},data] = xra%s(cellData,xraData,SOC,T);',upper(method));
      eval(execXRA);

      % Save tf names and x-location information
      if tt == 1, tfData = data; end
    end
  end

  % Convert cell array to a two dimensional structure
  % #T times #soc struct array with fields {xRA,T,SOC,A,B,C,D}
  out = reshape([out{:,:}],length(xraData.T),length(xraData.SOC));

  % Save pre-computed ROM models, cellData, xraData, tfData
  ROM = [];
  ROM.ROMmdls  = out;
  ROM.cellData = cellData;
  ROM.xraData  = xraData;
  ROM.tfData   = tfData;
end