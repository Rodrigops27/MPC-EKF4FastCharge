% function [ind, data] = loadParamTable(fileName, sheet)
% 
% Input:
%   fileName = full path to Excel spreadsheet of cell parameter values
%   sheet    = name of sheet within that spreadsheet to load
% Output:
%   ind      = structure containing indices (row numbers) of major sections
%   data     = data loaded from the sheet
%
% This utility function loads a sheet from an Excel spreadsheet and
% identifies the row positions of the main segments of a cell parameter
% table.
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

function [ind, data] = loadParamTable(fileName, sheet)
    
  % Read data from Excel sheet, convert to cell array.
  try
    % 'Format','auto' converts everything to text (much faster), only
    % needed for MATLAB >= R2020a.
    data = readtable(fileName,'Sheet',sheet,'ReadVariableNames',false,'Format','auto');
  catch ME
    if ~strcmp(ME.identifier,'MATLAB:table:parseArgs:BadParamName')
      rethrow(ME);
    end
    % 'Format' not an option in older versions of MATLAB; default
    % behavior OK.
    data = readtable(fileName,'Sheet',sheet,'ReadVariableNames',false);
  end
  data = table2cell(data);
  
  % Convert strings to numbers or NaN (emulate output of xlsread for
  % compatability).
  data(cellfun(@(x)all(ismissing(x)),data)) = {NaN};
  numeric = cellfun(@(x)str2double(x),data);
  indNumeric = ~isnan(numeric);
  data(indNumeric) = num2cell(numeric(indNumeric));

  % Determine row ranges for the main segments of the XLSX data file
  % (they don't need to be in any specific order)
  col1 = data(:,1); 
  cmpGen = strcmpi(col1,'#general'); indGen = find(cmpGen == 1,1);
  cmpCst = strcmpi(col1,'#const');   indCst = find(cmpCst == 1,1);
  cmpNeg = strcmpi(col1,'#neg');     indNeg = find(cmpNeg == 1,1);
  cmpSep = strcmpi(col1,'#sep');     indSep = find(cmpSep == 1,1);
  cmpPos = strcmpi(col1,'#pos');     indPos = find(cmpPos == 1,1);
  if isempty(indGen) || isempty(indCst) || isempty(indNeg) || ...
    isempty(indSep) || isempty(indPos)
    error(['Missing "#general", "#const", "#neg", "#sep", and/or ' ...
         '"#pos" section(s) in xslx file: '],fileName);
  end
  rowBounds = sort([indGen indCst indNeg indSep indPos length(col1)+1]);
  ind.gen = indGen:(rowBounds(find(rowBounds==indGen)+1)-1);
  ind.const = indCst:(rowBounds(find(rowBounds==indCst)+1)-1);
  ind.neg = indNeg:(rowBounds(find(rowBounds==indNeg)+1)-1);
  ind.sep = indSep:(rowBounds(find(rowBounds==indSep)+1)-1);
  ind.pos = indPos:(rowBounds(find(rowBounds==indPos)+1)-1);
end