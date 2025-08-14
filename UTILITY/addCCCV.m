% function FOM = addCCCV(FOM,Vmax)
% 
% Inputs:
%   FOM     = COMSOL object containing the full-order model, created by
%             genFOM.m
%   Vmax    = Maximum charge voltage, triggering CV mode
%
% Output:
%   FOM     = Updated COMSOL object containing the revised full-order model 
%
% This utility function modifies a COMSOL model to add a trigger that
% activates whenever the voltage exceeds Vmax; after that point in the
% simulation, the input current is modified to enforce constant-voltage
% output. This function uses the LiveLink for MATLAB interface. The
% returned COMSOL object can be saved to a file using mphsave.m, or loaded
% into the COMSOL GUI using mphlaunch.m. 
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

function FOM = addCCCV(FOM,Vmax)
  fprintf('Adding CC/CV trigger to model...\n');
  FOM.param.set('Vset', sprintf('%g',Vmax));
  FOM.component('mod1d').cpl.create('aveop1', 'Average'); % need
  FOM.component('mod1d').cpl('aveop1').selection.geom('geom1d', 0);
  FOM.component('mod1d').cpl('aveop1').selection.set(4);
  % FOM.physics('phi_s').feature('flux1').set('g', '-Iapp*xnorm');

  FOM.component('mod1d').physics.create('Iapp','GlobalEquations','geom1d');
  FOM.component('mod1d').physics('Iapp').identifier('Iapp');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set(...
    'DependentVariableQuantity', 'none');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set(...
    'CustomDependentVariableUnit', 'A');
  FOM.component('mod1d').physics.create('CCCV', 'Events', 'geom1d');
  FOM.component('mod1d').physics('CCCV').identifier('CCCV');
  FOM.component('mod1d').physics('CCCV').create('ds1','DiscreteStates',-1);
  FOM.component('mod1d').physics('CCCV').create('is1','IndicatorStates',-1);
  FOM.component('mod1d').physics('CCCV').create('impl1','ImplicitEvent',-1);

  FOM.component('mod1d').physics('Iapp').label('Compute CCCV');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set('name', 'Iapp');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set('equation', ...
    'CC*((Iapp-Ides)/Ides)+(1-CC)*(aveop1(phi_s)-Vset)/Vset');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set(...
    'initialValueU', 'Ides');
  FOM.component('mod1d').physics('CCCV').feature('ds1').set('dim', 'CC');
  FOM.component('mod1d').physics('CCCV').feature('ds1').set('dimInit', 1);
  FOM.component('mod1d').physics('CCCV').feature('ds1').set('dimDescr',...
    'Constant Current Flag');
  FOM.component('mod1d').physics('CCCV').feature('ds1').label(...
    'Initial States');
  FOM.component('mod1d').physics('CCCV').feature('is1').set('indDim','MaxV');
  FOM.component('mod1d').physics('CCCV').feature('is1').set('g',...
    'aveop1(phi_s)-Vset');
  FOM.component('mod1d').physics('CCCV').feature('is1').set('dimInit', 0);
  FOM.component('mod1d').physics('CCCV').feature('is1').set('dimDescr',...
    'Exceeded maximum voltage?');
  FOM.component('mod1d').physics('CCCV').feature('is1').label('Condition');
  FOM.component('mod1d').physics('CCCV').feature('impl1').set(...
    'condition','MaxV>0');
  FOM.component('mod1d').physics('CCCV').feature('impl1').set(...
    'reInitName', 'CC');
  FOM.component('mod1d').physics('CCCV').feature('impl1').set(...
    'reInitValue', 0);
  FOM.component('mod1d').physics('CCCV').feature('impl1').label(...
    'Reinitialize States');
