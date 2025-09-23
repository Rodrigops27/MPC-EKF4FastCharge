function [Vcell, obs, cellState] = OB_step(Iapp, Tc, cellState, ROM, initCfg)
% OB_STEP  One-timestep Output-Blend ROM simulation for MPC loops
%
%   [Vcell, obs, cellState] = OB_step(Iapp, Tc, cellState, ROM, initCfg)
%
% Inputs
%   Iapp      : applied current at this step (A)  (charging < 0)
%   Tc        : cell temperature in degC at this step
%   cellState : persistent state structure from previous call (pass [] on first call)
%   ROM       : xRA ROM structure (same used by outBlend)
%   initCfg   : struct used on first call only with fields:
%                 .SOC0   initial SOC in % (e.g., 10)
%                 .warnOff (optional) if true, suppresses short warnings
%
% Outputs
%   Vcell     : cell voltage at this step (V)
%   obs       : struct with per-step outputs (Ifdl/If/Idl/Phis/Phise/Phie/Thetae,
%               thetass, SOCn/ SOCp / cellSOC)
%   cellState : updated persistent state to feed next OB_step call
%
% Notes
%   - This is a single-step version of "outBlend" with the same blending,
%     nonlinear corrections, and voltage assembly. It keeps ROM states and
%     electrode-average SOC internally and advances them each time you call it.
%   - Sampling time is taken from ROM.xraData.Tsamp.
%
% Reference implementation adapted from "outBlend.m". 
% (Functions mirrored: setupIndsLocs, setupBlend, simStep, shortWarn)
%   Prof. Gregory L. Plett and Prof. M. Scott Trimboli
%   University of Colorado Colorado Springs (UCCS)
% as part of the Physics-Based Reduced-Order Model framework for lithium-ion
% batteries (see: Battery Management Systems, Volume III: Physics-Based
% Methods, Artech House, 2024).
% License: This file is distributed under the Creative Commons Attribution-ShareAlike
%   4.0 International License (CC BY-SA 4.0).
% -------------------------------------------------------------------------

  % -- One-time setup / first call ----------------------------------------
  if nargin < 3 || isempty(cellState) || ~isfield(cellState,'initialized') || ~cellState.initialized
    if nargin < 5 || ~isfield(initCfg,'SOC0')
      error('First call requires initCfg.SOC0 (in %).');
    end
    cellState = struct();
    cellState.initialized = true;
    cellState.warnOff = (nargin>=5) && isfield(initCfg,'warnOff') && logical(initCfg.warnOff);

    % Cache ROM bits
    cellState.ROMmdls  = ROM.ROMmdls;   % precomputed state-space models
    cellState.cellData = ROM.cellData;  % cell parameters
    cellState.Ts       = ROM.xraData.Tsamp; % sampling time (s)

    % Build indices/locations and blending matrices
    cellState.ROM = ROM; % keep full ROM for index building
    setupIndsLocs;
    [bigA,bigX,Tspts,Zspts,ZZ] = setupBlend;
    cellState.bigA  = bigA;
    cellState.bigX  = bigX;
    cellState.Tspts = Tspts;
    cellState.Zspts = Zspts;
    cellState.ZZ    = ZZ;

    % Initial electrode-average SOC (convert from % to local-electrode SOCs)
    Tk1 = Tc + 273.15;
    SOC0n = cellState.cellData.function.neg.soc(initCfg.SOC0/100, Tk1);
    SOC0p = cellState.cellData.function.pos.soc(initCfg.SOC0/100, Tk1);
    cellState.SOCnAvg = SOC0n;
    cellState.SOCpAvg = SOC0p;

    % Keep initial (for certain nonlinear terms as in outBlend)
    cellState.SOC0n = SOC0n;
    cellState.SOC0p = SOC0p;
  end

  % -- Do the single simulation step --------------------------------------
  [Vcell, newCellState, obs] = simStep(Iapp, Tc+273.15, cellState);

  % -- Export updated state for next call ----------------------------------
  cellState.bigX     = newCellState.bigX;
  cellState.SOCnAvg  = newCellState.SOCnAvg;
  cellState.SOCpAvg  = newCellState.SOCpAvg;

  % =========================================================================
  % Internal helpers mirrored/adapted from outBlend.m (single-step versions)
  % =========================================================================

  function setupIndsLocs
    ROMloc = cellState.ROM; % full ROM in state
    tfName = ROMloc.tfData.names;
    tfLocs = ROMloc.tfData.xLoc;
    cellState.tfLocs = tfLocs;

    % Negative electrode
    ind.negIfdl    = find(strcmp(tfName,'negIfdl') == 1);
    ind.negIf      = find(strcmp(tfName,'negIf') == 1);
    ind.negIdl     = find(strcmp(tfName,'negIdl') == 1);
    ind.negPhis    = find(strcmp(tfName,'negPhis') == 1);
    ind.negPhise   = find(strcmp(tfName,'negPhise') == 1);
    ind.negThetass = find(strcmp(tfName,'negThetass') == 1);

    % Positive electrode
    ind.posIfdl    = find(strcmp(tfName,'posIfdl') == 1);
    ind.posIf      = find(strcmp(tfName,'posIf') == 1);
    ind.posIdl     = find(strcmp(tfName,'posIdl') == 1);
    ind.posPhis    = find(strcmp(tfName,'posPhis') == 1);
    ind.posPhise   = find(strcmp(tfName,'posPhise') == 1);
    ind.posThetass = find(strcmp(tfName,'posThetass') == 1);

    % Electrolyte potential across cell width
    ind.negPhie = find(strcmp(tfName,'negPhie') == 1);
    ind.sepPhie = find(strcmp(tfName,'sepPhie') == 1);
    ind.posPhie = find(strcmp(tfName,'posPhie') == 1);
    loc.negPhie = tfLocs(ind.negPhie);
    loc.sepPhie = tfLocs(ind.sepPhie);
    loc.posPhie = tfLocs(ind.posPhie);

    ind.Phie = [ind.negPhie; ind.sepPhie; ind.posPhie];
    loc.Phie = [loc.negPhie; loc.sepPhie; loc.posPhie];

    % Electrolyte normalized concentration across cell width
    ind.negThetae  = find(strcmp(tfName,'negThetae')== 1);
    ind.sepThetae  = find(strcmp(tfName,'sepThetae')== 1);
    ind.posThetae  = find(strcmp(tfName,'posThetae')== 1);
    loc.negThetaes = tfLocs(ind.negThetae);
    loc.sepThetaes = tfLocs(ind.sepThetae);
    loc.posThetaes = tfLocs(ind.posThetae);

    ind.Thetae = [ind.negThetae; ind.sepThetae; ind.posThetae];
    loc.Thetae = [loc.negThetaes; loc.sepThetaes; loc.posThetaes];

    % Current-collector specific indices
    ind.negIfdl0    = find(strcmp(tfName,'negIfdl') == 1 & tfLocs == 0);
    ind.posIfdl3    = find(strcmp(tfName,'posIfdl') == 1 & tfLocs == 3);
    ind.negIf0      = find(strcmp(tfName,'negIf')   == 1 & tfLocs == 0);
    ind.posIf3      = find(strcmp(tfName,'posIf')   == 1 & tfLocs == 3);
    ind.negThetass0 = find(strcmp(tfName,'negThetass') == 1 & tfLocs == 0);
    ind.posThetass3 = find(strcmp(tfName,'posThetass') == 1 & tfLocs == 3);
    ind.negPhise0   = find(strcmp(tfName,'negPhise')   == 1 & tfLocs == 0);

    % Checks (same as outBlend)
    if isempty(ind.negIfdl0), error('Requires ifdl at negative collector'); end
    if isempty(ind.posIfdl3), error('Requires ifdl at positive collector'); end
    if isempty(ind.negIf0),   error('Requires if at negative collector');   end
    if isempty(ind.posIf3),   error('Requires if at positive collector');   end
    if loc.Thetae(1) > 0, error('Requires thetae at negative collector'); end
    if or(loc.Thetae(end)>3+eps,loc.Thetae(end)<3-eps)
      error('Requires thetae at positive collector'); 
    end
    if isempty(ind.negThetass0), error('Requires thetass at negative collector'); end
    if isempty(ind.posThetass3), error('Requires thetass at positive collector'); end
    if isempty(ind.negPhise0),   error('Requires phise at negative collector');   end
    if loc.Phie(1) == 0
      shortWarn('First phie x-location should not be zero. Ignoring');
      ind.Phie = ind.Phie(2:end);
      loc.Phie = loc.Phie(2:end);
    end
    if or(loc.Phie(end)>3+eps,loc.Phie(end)<3-eps)
      error('Requires phie at positive collector');
    end

    % stash
    cellState.ind = ind;
    cellState.loc = loc;
  end

  function [bigA,bigX,Tspts,Zspts,ZZ] = setupBlend
    Tspts = sort(ROM.xraData.T)+273.15;
    Zspts = sort(ROM.xraData.SOC/100);
    [TT,ZZ] = size(cellState.ROMmdls);

    % Column-wise diag(A) stack for each (T,Z) setpoint
    bigA = zeros(length(cellState.ROMmdls(1,1).A), TT*ZZ);
    for tt = 1:TT
      for zz = 1:ZZ
        indcol  = (tt-1)*ZZ + zz;
        ROMi    = cellState.ROMmdls(tt,zz);
        bigA(:,indcol) = real(diag(ROMi.A));
        % remove res0 from Phise rows (same as outBlend)
        ROMi.C(cellState.ind.negPhise,end) = 0;
        if any(cellState.ind.posPhise)
          ROMi.C(cellState.ind.posPhise,end) = 0;
        end
        cellState.ROMmdls(tt,zz) = ROMi; % store back
      end
    end
    bigX = zeros(size(bigA));
  end

  function [Vcell,newCellState,obs] = simStep(Iapp,T,cs)
    % Aliases
    ind      = cs.ind;
    loc      = cs.loc;
    Ts       = cs.Ts;
    ROMmdls  = cs.ROMmdls;
    bigA     = cs.bigA;
    bigX     = cs.bigX;
    Tspts    = cs.Tspts;
    Zspts    = cs.Zspts;
    ZZ       = cs.ZZ;
    cellData = cs.cellData;
    SOC0n    = cs.SOC0n; % used for wDL/Cdl like original code
    SOC0p    = cs.SOC0p;

    F         = cellData.const.F;
    R         = cellData.const.R;
    Q         = cellData.function.const.Q();
    Rc        = cellData.function.const.Rc();
    theta0n   = cellData.function.neg.theta0();
    theta0p   = cellData.function.pos.theta0();
    theta100n = cellData.function.neg.theta100();
    theta100p = cellData.function.pos.theta100();

    wDLn      = cellData.function.neg.wDL(SOC0n,T);
    wDLp      = cellData.function.pos.wDL(SOC0p,T);
    Cdln      = cellData.function.neg.Cdl(SOC0n,T);
    Cdlp      = cellData.function.pos.Cdl(SOC0p,T);
    nDLn      = cellData.function.neg.nDL();
    nDLp      = cellData.function.pos.nDL();
    Cdleffn   = (Cdln^(2-nDLn))*(wDLn^(nDLn-1));
    Cdleffp   = (Cdlp^(2-nDLp))*(wDLp^(nDLp-1));

    % Load present model state
    SOCnAvg = cs.SOCnAvg;
    SOCpAvg = cs.SOCpAvg;

    % ---- Step 1: Update (electrode-average) SOC & compute cell SOC
    obs.negSOC  = SOCnAvg;
    obs.posSOC  = SOCpAvg;
    obs.cellSOC = (SOCnAvg - theta0n)/(theta100n - theta0n);

    % Compute integrator input gains
    dUocpnAvg = cellData.function.neg.dUocp(SOCnAvg,T);
    dUocppAvg = cellData.function.pos.dUocp(SOCpAvg,T);
    dQn       = abs(theta100n-theta0n);
    dQp       = abs(theta100p-theta0p);
    res0n     = -dQn/(3600*Q - Cdleffn*dQn*dUocpnAvg);
    res0p     =  dQp/(3600*Q - Cdleffp*dQp*dUocppAvg);

    % Compute average negative-electrode and positive-electrode SOCs
    SOCnAvg = SOCnAvg + res0n*Iapp*Ts;
    SOCpAvg = SOCpAvg + res0p*Iapp*Ts;
    if SOCnAvg < 0, shortWarn('Average SOCn < 0'); SOCnAvg = 0; end
    if SOCnAvg > 1, shortWarn('Average SOCn > 1'); SOCnAvg = 1; end
    if SOCpAvg < 0, shortWarn('Average SOCp < 0'); SOCpAvg = 0; end
    if SOCpAvg > 1, shortWarn('Average SOCp > 1'); SOCpAvg = 1; end

    % ---------------------------------------------------------------------
    % Step 2: Pick 4 nearest (T,SOC) setpoints; compute yk at those
    % Z upper/lower
    %         y[k] = C*x[k] + D*iapp[k]
    % ---------------------------------------------------------------------
    [~,iZ] = sort(abs(obs.cellSOC - Zspts));
    if numel(iZ) > 1
      iZupper = max(iZ(1), iZ(2));
      iZlower = min(iZ(1), iZ(2));
      Zupper  = Zspts(iZupper);
      Zlower  = Zspts(iZlower);
    else
      iZupper=iZ; iZlower=iZ; Zupper=Zspts; Zlower=Zspts;
    end
    % T upper/lower
    [~,iT] = sort(abs(T - Tspts));
    if numel(iT) > 1
      iTupper = max(iT(1), iT(2));
      iTlower = min(iT(1), iT(2));
      Tupper  = Tspts(iTupper);
      Tlower  = Tspts(iTlower);
    else
      iTupper=iT; iTlower=iT; Tupper=Tspts; Tlower=Tspts;
    end

    % Four corner outputs
    yk1 = ROMmdls(iTlower,iZlower).C*bigX(:,(iTlower-1)*ZZ+iZlower) + ROMmdls(iTlower,iZlower).D*Iapp;
    yk2 = ROMmdls(iTlower,iZupper).C*bigX(:,(iTlower-1)*ZZ+iZupper) + ROMmdls(iTlower,iZupper).D*Iapp;
    yk3 = ROMmdls(iTupper,iZlower).C*bigX(:,(iTupper-1)*ZZ+iZlower) + ROMmdls(iTupper,iZlower).D*Iapp;
    yk4 = ROMmdls(iTupper,iZupper).C*bigX(:,(iTupper-1)*ZZ+iZupper) + ROMmdls(iTupper,iZupper).D*Iapp;

    % ---- Step 3: Advance states for *all* precomputed ROMs (A diagonal, B=1)
    bigX = cs.bigA.*bigX + Iapp;

    % ---- Step 4: Blend linearly to present setpoint
    alphaZ = 0; alphaT = 0;
    if Zupper ~= Zlower, alphaZ = (obs.cellSOC - Zlower)/(Zupper - Zlower); end
    if Tupper ~= Tlower, alphaT = (T - Tlower)/(Tupper - Tlower); end
    
    yk = (1-alphaT)*((1-alphaZ)*yk1 + alphaZ*yk2) + alphaT*((1-alphaZ)*yk3 + alphaZ*yk4);

    % ---- Step 5: Nonlinear corrections + assemble observables
    % ifdl
    obs.negIfdl   = yk(ind.negIfdl);
    obs.posIfdl   = yk(ind.posIfdl);
    obs.negIfdl0  = yk(ind.negIfdl0);
    obs.posIfdl3  = yk(ind.posIfdl3);
    % if
    obs.negIf     = yk(ind.negIf);
    obs.posIf     = yk(ind.posIf);
    obs.negIf0    = yk(ind.negIf0);
    obs.posIf3    = yk(ind.posIf3);
    % idl
    obs.negIdl    = yk(ind.negIdl);
    obs.posIdl    = yk(ind.posIdl);
    % thetass (add average SOC offsets)
    obs.negThetass = yk(ind.negThetass) + SOC0n;
    obs.posThetass = yk(ind.posThetass) + SOC0p;
    obs.negThetass0 = yk(ind.negThetass0) + SOC0n;
    obs.posThetass3 = yk(ind.posThetass3) + SOC0p;
    % clamp
    obs.negThetass  = min(max(obs.negThetass,  1e-6), 1-1e-6);
    obs.posThetass  = min(max(obs.posThetass,  1e-6), 1-1e-6);
    obs.negThetass0 = min(max(obs.negThetass0, 1e-6), 1-1e-6);
    obs.posThetass3 = min(max(obs.posThetass3, 1e-6), 1-1e-6);

    % phise (add Uocp(avg))
    UocpnAvg = cellData.function.neg.Uocp(obs.negSOC, T);
    UocppAvg = cellData.function.pos.Uocp(obs.posSOC, T);
    obs.negPhise  = yk(ind.negPhise) + UocpnAvg;
    obs.posPhise  = yk(ind.posPhise) + UocppAvg;
    obs.negPhise0 = yk(ind.negPhise0) + UocpnAvg;

    % Phie: first at x=0, then across width
    obs.Phie = zeros(numel(ind.Phie)+1,1);
    obs.Phie(1)      = -obs.negPhise0;
    obs.Phie(2:end)  = yk(ind.Phie) - obs.negPhise0;

    % thetae
    obs.Thetae = yk(ind.Thetae) + 1;
    obs.Thetae = max(obs.Thetae, 1e-6);

    % Overpotentials (asinh method)
    k0n = cellData.function.neg.k0(obs.negSOC, T);
    k0p = cellData.function.pos.k0(obs.posSOC, T);
    i0n = k0n*sqrt(obs.Thetae(1)*(1-obs.negThetass0)*obs.negThetass0);
    i0p = k0p*sqrt(obs.Thetae(end)*(1-obs.posThetass3)*obs.posThetass3);
    negEta0 = 2*R*T/F*asinh(obs.negIf0/(2*i0n));
    posEta3 = 2*R*T/F*asinh(obs.posIf3/(2*i0p));

    % Cell voltage
    Uocpn0 = cellData.function.neg.Uocp(obs.negThetass0, T);
    Uocpp3 = cellData.function.pos.Uocp(obs.posThetass3, T);
    Rfn    = cellData.function.neg.Rf(obs.negSOC, T);
    Rfp    = cellData.function.pos.Rf(obs.posSOC, T);
    Vcell  = posEta3 - negEta0 + yk(ind.Phie(end)) + Uocpp3 - Uocpn0 ...
             + (Rfp*obs.posIfdl3 - Rfn*obs.negIfdl0);
    % ohmic Rc
    Vcell  = Vcell - Rc*Iapp;

    % Solid potentials
    obs.negPhis = yk(ind.negPhis);
    obs.posPhis = yk(ind.posPhis) + Vcell;

    % Return updated ROM states & electrode averages
    newCellState.bigX    = bigX;
    newCellState.SOCnAvg = SOCnAvg;
    newCellState.SOCpAvg = SOCpAvg;

    % Also helpful to return commonly-used scalars
    obs.Vcell = Vcell;
  end

  function shortWarn(msg)
    if cellState.warnOff, return; end
    persistent warnState
    if strcmpi(msg,'on')
      warnState = [];
    elseif strcmpi(msg,'off')
      warnState = 1;
    elseif isempty(warnState)
      % Minimal, line-number-free warning like the original
      fprintf(2, ' - Warning: %s\n', msg);
    end
  end
end
