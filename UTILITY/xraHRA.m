% function [ROM,tfData] = xraHRA(cellData,xraData,SOC,T)
%
% Inputs:
%   cellData = cell-parameter data structure loaded with "loadCellParams"
%   xraData  = xra-control parameters loaded with "loadXRA"
%   SOC      = cell state of charge setpoint (between 0 and 1) for the
%              model to be created
%   T        = cell temperature (in K) for the model to be created
%
% Outputs:
%   ROM      = a structure containing the ROM data
%
% This function creates a reduced-order model for the cell defined by
% "cellData" using the XRA tuning parameters in "xraData" for setpoint
% defined by (SOC,T)... using the hybrid realization algorithm (HRA).
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

function [ROM,tfData] = xraHRA(cellData,xraData,SOC,T)
% ----------------- Load HRA tuning parameters -----------------
Ts   = xraData.Tsamp;    % sampling period [s]
n1   = xraData.n1;       % initial # non-integrator poles
n    = xraData.n;        % final # non-integrator poles
Gc   = xraData.HRA_Gc;   % alias cutoff magnitude (dB)
ii   = xraData.HRA_ii;   % self-adjusts...require n <= (ii+1)*numOut
tflist = xraData.tf;

wnyquist = pi/Ts;          % the Nyquist sampling frequency [rad/s]
W = unique(xraData.HRA_W)*wnyquist; % up to the Nqyuist rate
npts = length(W);          % Need to be >> n

% WH: Define number of inputs (single iapp for LIB cells).
numIn = 1;

% ----------------- Find discrete-time freq. resp. -----------------
% WH: automatically determine number of aliases needed for each TF.
cellData0 = evalSetpoint(cellData,0,SOC,T);
[kmaxlist,wclist,tflistFlat] = autoKmax(tflist, wnyquist, cellData0, Gc);

% WH: ensure at least HRA_Kmax aliases are used.
Kmax_min = xraData.HRA_Kmax_min;
kmaxlist(kmaxlist<Kmax_min) = Kmax_min;
if isfield(xraData,'debug') && xraData.debug
  wclist = sprintf('%.3g ',wclist);
  kmaxstr = sprintf('%d ',kmaxlist);
  fprintf(' - wc  [r/s]: %s\n',wclist);
  fprintf(' - Kmax     : %s\n',kmaxstr);
end

numOut = length(tflistFlat);
numTFs = numOut;  % number of TFs before removing zero TFs below

% Check if ii is large enough
if (ii+1)*numOut < n1
  warning('ii not large enough... adjusting');
  ii = ceil(n1/numOut - 1);
end

% WH: perform C2D for each TF.
Gdest = [];    % DT frequency response (dim1=output, dim2=frequency)
hfGain = [];   % D term vector
res0 = [];     % integrator residuals vector
K = [];        % dc gain vector
tfData = struct('names',[],'xLoc',[]);
isZero = false(numOut,1);
for ktf = 1:numOut
  tf = tflistFlat(ktf);
  kmax = kmaxlist(ktf);

  % Compute continuous-time frequencies needed for C2D.
  M = (ones(2*kmax+1,1)*W) + 2*pi/Ts*((-kmax:kmax)'*ones(1,length(W)));
  s = 1j*M(:).';

  % Evaluate model parameters at the setpoint
  cellData_ = evalSetpoint(cellData,s,SOC,T);

  % Evaluate frequency response of the current TF.
  [freqResp_,hfGain_,res0_,tfData_] = evalTF(tf,s,cellData_);

  % Remove D term for now; add back into final model.
  freqResp_ = (freqResp_ - hfGain_*ones(1,(2*kmax+1)*npts)).';

  % If this TF is identically zero, delete it for now.
  zeroTest = sum(abs(freqResp_));
  isZero(ktf) = zeroTest == 0;
  if isZero(ktf)
    continue;  % skip rest of this loop
  end

  % Set dc gain to 1 (add back into final model)
  dcGain = real(freqResp_(M(:)==0,:)); % This is dcgain without D term
  dcGain(dcGain==0) = 1;
  dcGain(isinf(dcGain)) = 1;
  dcGain(isnan(dcGain)) = 1;
  freqResp_(M(:)==0,:) = dcGain;
  K_ = dcGain;
  freqResp_ = freqResp_/K_;

  % Add together frequency alises to create discrete-time freq resp
  newData = freqResp_(:).*(1-exp(-1j*M(:)*Ts))./(1j*M(:));
  newData(M(:)==0) = freqResp_(M(:)==0)*Ts;
  geval = reshape(newData.',size(M));
  gdest = (1/Ts)*sum(geval);
  gdest = reshape(gdest.',length(W),1).';

  % Store results.
  Gdest = [Gdest; gdest]; %#ok<AGROW>
  hfGain = [hfGain; hfGain_]; %#ok<AGROW>
  res0 = [res0; res0_]; %#ok<AGROW>
  tfData.names = [tfData.names; tfData_.names];
  tfData.xLoc = [tfData.xLoc; tfData_.xLoc];
  K = [K; K_]; %#ok<AGROW>

  % WH: debug continuous-to-discrete frequency response conversion.
  if isfield(xraData,'debugC2D') && xraData.debugC2D
    M__ = (ones(kmax+1,1)*W) + 2*pi/Ts*((0:kmax)'*ones(1,length(W)));
    W__ = sort(M__(:).');
    [freqResp__,hfGain_] = evalTF(tf,1j*W__.',cellData_);
    freqResp__ = (freqResp__ - hfGain_*ones(1,(kmax+1)*npts)).';
    dcGain__ = real(freqResp__(W__==0,:)); % This is dcgain without D term
    dcGain__(dcGain__==0) = 1;
    dcGain__(isinf(dcGain__)) = 1;
    dcGain__(isnan(dcGain__)) = 1;
    dcGain__(W__==0,:) = dcGain__;
    K__ = dcGain__;
    freqResp__ = freqResp__/K__;
    % Plotting.
    figure;
    subplot(1,3,1);
    semilogx(W__,20*log10(abs(freqResp__(:))),'b'); hold on;
    semilogx(W,20*log10(abs(gdest(:))),'r.');
    xline(wclist(ktf),'m');
    title('Bode Magnitude');
    xlabel('Radian frequency [rad/s]');
    ylabel('20log|G(\omega)|');
    legend('Continuous','Discrete','Location','best');
    subplot(1,3,2);
    semilogx(W__,angle(freqResp__(:))*180/pi,'b'); hold on;
    semilogx(W,angle(gdest(:))*180/pi,'r.');
    xlabel('Radian frequency [rad/s]');
    ylabel('\angle G(\omega) [deg]');
    title('Bode Phase');
    subplot(1,3,3);
    plot(real(freqResp__(:)),imag(freqResp__(:)),'b'); hold on;
    plot(real(gdest(:)),imag(gdest(:)),'r.'); hold on;
    xlabel('G''(\omega)');
    ylabel('G''''(\omega)');
    title('Nyquist');
    setAxesNyquist;
    tfID = sprintf( ...
      'TF # = %d, TF name = %s, loc = %d', ...
      ktf,tfData_.names{1},tfData_.xLoc(1) ...
      );
    plttitle = sprintf( ...
      '%s: continuous-to-discrete conversion (kmax=%d)',tfID,kmax ...
      );
    sgtitle(plttitle);
    if exist('thesisFormat', 'file')
      try
        thesisFormat('PlotBoxPaddingInches', [0 0 0 0.2]);
      catch
        thesisFormat([0 0 0 0.2]);
      end % try
    end % if
  end
end % for

% Convert dc gain vector to diagonal matrix.
K = diag(K);

% Update number of outputs to account for deleted TFs.
indZero = find(isZero); %#ok<EFIND>
indNonzero = find(~isZero);
numOut = size(Gdest,1);

% WH: debug continuous-to-discrete frequency response conversion with
% discrete-time KK test.
if isfield(xraData,'debugC2DKK') && xraData.debugC2DKK
  Omega_ = W*Ts;
  N_ = npts/40;
  for kk = 1:numOut
    G_ = Gdest(kk,:);
    kkData = dtKK(Omega_(Omega_~=0),G_(Omega_~=0),'minM',N_,'maxM',N_);
    figure;
    subplot(1,2,1);
    plot(real(G_),imag(G_),'b.'); hold on;
    plot(real(kkData.Gmodel),imag(kkData.Gmodel),'g.');
    xlabel("G'(e^{j\Omega})");
    ylabel("G''(e^{j\Omega})");
    legend('Truth','KK Model','Location','best');
    title('Nyquist: KK Test');
    setAxesNyquist;
    subplot(1,2,2);
    semilogx(Omega_(Omega_~=0),abs(kkData.residuals)*100);
    xlabel('\Omega [rad/sample]');
    ylabel('KK Residual [%]');
    title('Residual: KK Test');
    tfID = sprintf( ...
      'TF # = %d, TF name = %s, loc = %d', ...
      kk,tfData.names{kk},tfData.xLoc(kk) ...
      );
    plttitle = sprintf( ...
      '%s: KK test on DT frequency response (kmax=%d)',tfID,kmaxlist(kk) ...
      );
    sgtitle(plttitle);
    if exist('thesisFormat', 'file')
      try
        thesisFormat('PlotBoxPaddingInches', [0 0 0 0.3]);
      catch
        thesisFormat([0 0 0 0.3]);
      end % try
    end % if
  end % for numOut
end % if debug

% ------------------------ Execute HRA ------------------------
% At this point, W is the frequency vector; Gdest is Gd(exp(j*w*Ts)).
% Need to assemble Y = [Gd; exp(j*w*Ts)*Gd; exp(2j*w*Ts)*Gd; ...]
% and U = [exp(0j*w*Ts); exp(1j*w*Ts); exp(2j*w*Ts); ...]
% Then, find [U;G] = [L11 0; L21 L22] * [Q1; Q2]
% Then, [U1 U2] * [S11 0; 0 0] * [V1'; V2'] = L22
% Then, Ok = U1*sqrt(S11)
% Then, C = Ok(1:p,1:n)
% Finally, [Ok(1:p*(k-1),1:n)] * A = [Ok(p+1:k*p,1:n)].

% Make data matrices
% WH: support multiple inputs (numIn>1)
Ucell = cell(ii+1, npts);
Ycell = cell(ii+1, npts);
for ki = 1:ii+1
  ind = ki - 1;  % 0-based oversizing index
  for kw = 1:npts
    w0 = W(kw);        % radian frequency
    G0 = Gdest(:,kw);  % frequency response at w0
    Ucell{ki, kw} = exp(1j*w0*ind*Ts)*eye(numIn);
    Ycell{ki, kw} = exp(1j*w0*ind*Ts)*G0;
  end % for w
end % for i
U = cell2mat(Ucell);
Y = cell2mat(Ycell);

% Perform LQ decomposition of [U;Y]
% WH: support multiple inputs (numIn>1)
Ur = [real(U) imag(U)];
Yr = [real(Y) imag(Y)];
LQmat = [Ur; Yr];
[Q,R2] = qr(LQmat',"econ");
L2 = R2';
Q = Q';

% Partition the L matrix.
% WH: correct partitioning.
[nru, ncu] = size(U); % number of rows/cols in U (nc(U)=nc(Ur)/2)
if size(LQmat,1) < size(LQmat,2)
  % We have the economy-size LQ decomposition.
  rpi = nru;  % row partition index
  cpi = nru;  % col partition index
else
  % We have the full-size LQ decomposition.
  rpi = nru;  % row partition index
  cpi = ncu;  % col partition index
end
L22 = L2(rpi+1:end,cpi+1:end);

% WH: Ensure the LQ decomposition is doing what we think it's doing.
if isfield(xraData,'debug') && xraData.debug
  L11 = L2(1:rpi,1:cpi);
  L12_zero = L2(1:rpi,cpi+1:end);  % should be 0 matrix
  L21 = L2(rpi+1:end,1:cpi);
  Q1 = Q(1:cpi,:);
  Q2 = Q(cpi+1:end,:);

  % Check dimensional conformality.
  assert( ...
    size(L11,2)==size(Q1,1), ...
    "Cannot multiply L11 (%dx%d) and Q1 (%dx%d).", ...
    size(L11,1), size(L11,2), size(Q1,1), size(Q1,2) ...
    );
  assert( ...
    all(size(L11*Q1)==size(Ur)), ...
    "L11*Q1 (%dx%d) is not dimensionally conformable with Ur (%dx%d).", ...
    size(L11*Q1,1), size(L11*Q1,2), size(Ur,1), size(Ur,2) ...
    );
  assert( ...
    size(L21,2)==size(Q1,1), ...
    "Cannot multiply L21 (%dx%d) and Q1 (%dx%d).", ...
    size(L21,1), size(L21,2), size(Q1,1), size(Q1,2) ...
    );
  assert( ...
    all(size(L21*Q1)==size(Yr)), ...
    "L21*Q1 (%dx%d) is not dimensionally conformable with Yr (%dx%d).", ...
    size(L21*Q1,1), size(L21*Q1,2), size(Yr,1), size(Yr,2) ...
    );
  assert( ...
    size(L22,2)==size(Q2,1), ...
    "Cannot multiply L22 (%dx%d) and Q2 (%dx%d).", ...
    size(L22,1), size(L22,2), size(Q2,1), size(Q2,2) ...
    );
  assert( ...
    all(size(L22*Q2)==size(Yr)), ...
    "L22*Q2 is not dimensionally conformable with Yr.", ...
    size(L22*Q2,1), size(L22*Q2,2), size(Yr,1), size(Yr,2) ...
    );

  % Check upper right block of L is 0.
  assert( ...
    all(L12_zero==0,"all"), ...
    "Partitioned L matrix does not have 0 as L12." ...
    );

  % Check that Q1 and Q2 are orthonormal.
  assert( ...
    norm(Q1*Q1'-eye(size(Q1*Q1')),"fro")/ncu<=1e-10, ...
    "Q1 is not orthonormal." ...
    );
  assert( ...
    norm(Q2*Q2'-eye(size(Q2*Q2')),"fro")/ncu<=1e-10, ...
    "Q2 is not orthonormal." ...
    );
end

% Find L22 via the LQ decomposition
[U,S,~] = svd(L22);
U1 = U(1:end,1:n1); S11 = S(1:n1,1:n1);
Ok = real(U1*sqrt(S11));

% Form A (order n),and D matrices
A = Ok(1:numOut*ii,1:n1)\Ok(numOut+1:(ii+1)*numOut,1:n1);
D = hfGain;

% Diagonalize A and replace unstable poles by their reciprocals
% and complex poles by their magnitude
eigA = eig(A);
eigA(abs(eigA)>1) = 1./eigA(abs(eigA)>1); % flip unstable
eigA = unique(sort(eigA,'descend'));
A = diag(abs(eigA)); % make sure real and positive
diagA = diag(A);

% ----------------- Compute C matrix (order n) -----------------
B = ones(length(diagA),1); % order n

% Use Lagrange multiplier method to compute C to enforce dc gain
C = zeros(numOut,length(diagA));
for kk = 1:numOut
  G0 = real(Gdest(kk,find(W==0,1,'first'))); % should be real...
  M1 = zeros(1,length(diagA));
  M2 = 0*A;
  for k = 1:length(W) % W,Gdest
    Mx = 1./(exp(1j*W(k)*Ts)-diagA);
    M1 = M1 + real(Gdest(kk,k)*Mx');
    M2 = M2 + real(Mx*Mx');
  end
  M0 = 1./(1-diagA);
  M  = [2*M2' -M0; M0' 0];
  N  = [2*M1'; G0'];
  CL = pinv(M)*N;
  C(kk,:) = real(CL(1:end-1,:))'; % should be real, but sometimes isn't

  % WH: debug C matrix optimization
  if (...
      isfield(xraData,'debugC') && xraData.debugC && ...
      isfield(xraData,'debugCfullOrder') && xraData.debugCfullOrder ...
      )
    H = zeros(size(W));
    for k = 1:length(W)
      Mx = diag(1./(exp(1j*W(k)*Ts)-diagA))*B;
      H(k) = C(kk,:)*Mx;
    end % for
    wpole = -log(diagA)/Ts;
    G = Gdest(kk,:);
    figure;
    subplot(1,3,1);
    semilogx(W,20*log10(abs(G)),'b'); hold on;
    semilogx(W,20*log10(abs(H)),'rx');
    xline(wpole,'m--');
    xlim([min(W(W~=0)) max(W)]);
    title('Bode Magnitude');
    xlabel('Radian frequency [rad/s]');
    ylabel('20log|G(\omega)|');
    legend('Truth','Estimate','Pole Frequencies','Location','best');
    subplot(1,3,2);
    semilogx(W,angle(G)*180/pi,'b'); hold on;
    semilogx(W,angle(H)*180/pi,'rx');
    xlabel('Radian frequency [rad/s]');
    ylabel('\angle G(\omega) [deg]');
    title('Bode Phase');
    subplot(1,3,3);
    plot(real(G),imag(G),'b'); hold on;
    plot(real(H),imag(H),'rx'); hold on;
    xlabel('G''(\omega)');
    ylabel('G''''(\omega)');
    title('Nyquist');
    setAxesNyquist;
    sgtitle('C matrix optimization (Lagrange multiplier method, full order)');
    if exist('thesisFormat', 'file')
      try
        thesisFormat('PlotBoxPaddingInches', [0 0 0 0.2]);
      catch
        thesisFormat([0 0 0 0.2]);
      end % try
    end % if
  end % if debugC
end % for

% WH: debug: compute C matrix with fmincon.
if (...
    isfield(xraData,'debugCfmin') && xraData.debugCfmin && ...
    isfield(xraData,'debugCfullOrder') && xraData.debugCfullOrder ...
    )
  Cfmin = ones(numOut,length(diagA));
  for kk = 1:numOut
    G0 = real(Gdest(kk,find(W==0,1,'first')));
    M0 = diag(1./(1-diagA))*B;
    G = Gdest(kk,:);
    Aeq = M0';
    beq = G0;
    [cost, opts] = getCmatrixCostFunction(A,B,G,W,Ts);
    Cfmin(kk,:) = fmincon(cost,Cfmin(kk,:),[],[],Aeq,beq,[],[],[],opts);
    % Plotting.
    wpole = -log(diagA)/Ts;
    [~, H] = cost(Cfmin(kk,:));
    figure;
    subplot(1,3,1);
    semilogx(W,20*log10(abs(G)),'b'); hold on;
    semilogx(W,20*log10(abs(H)),'rx');
    xline(wpole,'m--');
    xlim([min(W(W~=0)) max(W)]);
    title('Bode Magnitude');
    xlabel('Radian frequency [rad/s]');
    ylabel('20log|G(\omega)|');
    legend('Truth','Estimate','Pole Frequencies','Location','best');
    subplot(1,3,2);
    semilogx(W,angle(G)*180/pi,'b'); hold on;
    semilogx(W,angle(H)*180/pi,'rx');
    xlabel('Radian frequency [rad/s]');
    ylabel('\angle G(\omega) [deg]');
    title('Bode Phase');
    subplot(1,3,3);
    plot(real(G),imag(G),'b'); hold on;
    plot(real(H),imag(H),'rx'); hold on;
    xlabel('G''(\omega)');
    ylabel('G''''(\omega)');
    title('Nyquist');
    setAxesNyquist;
    sgtitle('C matrix optimization (fmincon method, full order)');
    if exist('thesisFormat', 'file')
      try
        thesisFormat('PlotBoxPaddingInches', [0 0 0 0.2]);
      catch
        thesisFormat([0 0 0 0.2]);
      end % try
    end % if
  end % for
end % if debugC

% ----------- Balanced reduction, eliminate complex poles; rescale C -------
sysHRA = ss(A,B,K*C,D,Ts);
opt = balredOptions('StateElimMethod','Truncate');
sysHRA = balred(sysHRA,n,opt);
sysHRA = canon(sysHRA,'modal',Inf);
[A,~,~,~] = ssdata(sysHRA);

% Again, diagonalize A and replace unstable poles by their reciprocals
% and complex poles by their magnitude
eigA = eig(A);
if isfield(xraData,'debug')
  if xraData.debug
    if any(abs(eigA) >= 1), cprintf([1,1/2,0],' - Unstable A... corrected\n'); end
    if any(real(eigA) < 0), cprintf([1,1/2,0],' - Oscillatory A... corrected\n'); end
    if any(eigA ~= conj(eigA)), cprintf([1,1/2,0],' - Complex A... corrected\n'); end
  end
end

eigA(abs(eigA)>1) = 1./eigA(abs(eigA)>1); % flip unstable
eigA = abs(sort(eigA,'descend'));
A = diag(eigA); % make sure real and positive
diagA = diag(A);

% WH: print pole break frequencies for debug.
if isfield(xraData,'debug') && xraData.debug
  strPoles = sprintf('%.3f ', diagA);
  fprintf(' - poles: %s\n', strPoles);
end

% ----------------- Compute C matrix again (order nfinal) ---------------
B = ones(n,numIn);  % WH: support multiple inputs

% Use Lagrange multiplier method to compute C to enforce dc gain
C = zeros(numOut,n);
for kk = 1:numOut
  G0 = real(Gdest(kk,find(W==0,1,'first'))); % should be real...
  M1 = zeros(1,n);
  M2 = 0*A;
  for k = 1:length(W) % W,Gdest
    Mx = 1./(exp(1j*W(k)*Ts)-diagA);
    M1 = M1 + real(Gdest(kk,k)*Mx');
    M2 = M2 + real(Mx*Mx');
  end
  M0 = 1./(1-diagA);
  M = [2*M2' -M0; M0' 0];
  N = [2*M1'; G0'];
  CL = pinv(M)*N;
  C(kk,:) = real(CL(1:end-1,:))'; % should be real, but sometimes isn't

  % WH: debug C matrix optimization
  if isfield(xraData,'debugC') && xraData.debugC
    G = Gdest(kk,:);
    H = zeros(size(W));
    for k = 1:length(W)
      Mx = diag(1./(exp(1j*W(k)*Ts)-diagA))*B;
      H(k) = C(kk,:)*Mx;
    end % for
    wpole = -log(diagA)/Ts;
    zz = sort(zero(ss(A,B,C(kk,:),0,Ts)),'descend');
    wzero = -log(abs(zz))/Ts;
    fprintf(' - (TF #%d) zeros: %s\n',kk,sprintf('%.3f ',zz));

    figure;
    subplot(1,3,1);
    semilogx(W,20*log10(abs(G)),'b'); hold on;
    semilogx(W,20*log10(abs(H)),'rx');
    xline(NaN,'m-'); xline(NaN,'k--');
    xline(wpole,'m-');
    xline(wzero,'k--');
    xlim([min(W(W~=0)) max(W)]);
    title('Bode Magnitude');
    xlabel('Radian frequency [rad/s]');
    ylabel('20log|G(e^{j\omega t_s})|');
    legend('Truth','Estimate','\omega_{pole}','\omega_{zero}','Location','best');
    subplot(1,3,2);
    semilogx(W,unwrap(angle(G))*180/pi,'b'); hold on;
    semilogx(W,unwrap(angle(H))*180/pi,'rx');
    xlabel('Radian frequency [rad/s]');
    ylabel('\angle G(e^{j\omega t_s}) [deg]');
    title('Bode Phase');
    subplot(1,3,3);
    plot(real(G),imag(G),'b'); hold on;
    plot(real(H),imag(H),'rx'); hold on;
    xlabel('G''(e^{j\omega t_s})');
    ylabel('G''''(e^{j\omega t_s})');
    title('Nyquist');
    setAxesNyquist;
    tfID = sprintf( ...
      'TF # = %d, TF name = %s, loc = %d', ...
      kk,tfData.names{kk},tfData.xLoc(kk) ...
      );
    plttitle = sprintf( ...
      ['%s: C matrix optimization ' ...
      '(Lagrange multiplier method, final order)'],tfID ...
      );
    sgtitle(plttitle);
    if exist('thesisFormat', 'file')
      try
        thesisFormat('PlotBoxPaddingInches', [0 0 0 0.2]);
      catch
        thesisFormat([0 0 0 0.2]);
      end % try
    end % if
  end % if debugC
end

% WH: debug: compute C matrix with fmincon.
if isfield(xraData,'debugCfmin') && xraData.debugCfmin
  Cfmin = ones(numOut,n);
  for kk = 1:numOut
    G0 = real(Gdest(kk,find(W==0,1,'first')));
    M0 = diag(1./(1-diagA))*B;
    G = Gdest(kk,:);
    Aeq = M0';
    beq = G0;
    [cost, opts] = getCmatrixCostFunction(A,B,G,W,Ts);
    Cfmin(kk,:) = fmincon(cost,Cfmin(kk,:),[],[],Aeq,beq,[],[],[],opts);
    % Plotting.
    wpole = -log(diagA)/Ts;
    [~, H] = cost(Cfmin(kk,:));
    figure;
    subplot(1,3,1);
    semilogx(W,20*log10(abs(G)),'b'); hold on;
    semilogx(W,20*log10(abs(H)),'rx');
    xline(wpole,'m--');
    xlim([min(W(W~=0)) max(W)]);
    title('Bode Magnitude');
    xlabel('Radian frequency [rad/s]');
    ylabel('20log|G(\omega)|');
    legend('Truth','Estimate','Pole Frequencies','Location','best');
    subplot(1,3,2);
    semilogx(W,angle(G)*180/pi,'b'); hold on;
    semilogx(W,angle(H)*180/pi,'rx');
    xlabel('Radian frequency [rad/s]');
    ylabel('\angle G(\omega) [deg]');
    title('Bode Phase');
    subplot(1,3,3);
    plot(real(G),imag(G),'b'); hold on;
    plot(real(H),imag(H),'rx'); hold on;
    xlabel('G''(\omega)');
    ylabel('G''''(\omega)');
    title('Nyquist');
    setAxesNyquist;
    sgtitle('C matrix optimization (fmincon method, final order)');
    if exist('thesisFormat', 'file')
      try
        thesisFormat('PlotBoxPaddingInches', [0 0 0 0.2]);
      catch
        thesisFormat([0 0 0 0.2]);
      end % try
    end % if
  end % for
end % if debugC

% Add back in the dc gain.
C = K*C;

% Add back in integrator ...
if ~all(res0 == 0)       % Note: res0 is cts-time residue, not disc-time
  A(n+1,n+1) = 1; B = [B;1]; C = [C res0*Ts]; %20211102: Add factor of Ts
end

% Add back in deleted TFs...
if ~isempty(indZero)
  Csmall = C; Dsmall = D;
  C = zeros(numTFs,size(Csmall,2));
  D = zeros(numTFs,1);
  C(indNonzero,:) = Csmall;
  D(indNonzero,:) = Dsmall;
end

% Save to structure
ROM.xRA = 'HRA';
ROM.T   = T;
ROM.SOC = SOC;
ROM.A   = A;
ROM.B   = B;
ROM.C   = C;
ROM.D   = D;

% WH: return analog radian frequency
ROM.omega = W;

% WH: return original discrete-time frequency response
ROM.Gdest = K*Gdest + diag(D)*ones(numOut,npts);

if isfield(xraData,'debug') && xraData.debug
  % plotting code... make system w/o integrator
  % WH: only remove integrator if present
  if ~all(res0 == 0)
    sysHRA = ss(A(1:end-1,1:end-1),B(1:end-1),C(:,1:end-1),D,Ts);
  else
    sysHRA = ss(A,B,C,D,Ts);
  end

  % Get frequency response of xra-realized system
  [magEst,phaseEst] = bode(sysHRA,W);
  magEst = squeeze(magEst); phaseEst = squeeze(phaseEst);
  % WH: special case when we have a single output.
  if numOut == 1
    magEst = magEst.';
    phaseEst = phaseEst.';
  end

  % Add dcgain and D term back into true freq response for plotting
  Gdest = K*Gdest + diag(D)*ones(numOut,npts);
  magTrue = abs(Gdest); phaseTrue = angle(Gdest);

  % Condition phaseEst, phaseTrue so they start in same multiple of 2*pi
  % Note that "unwrap" works on radian angles, hence the conversions...
  phaseTrue = unwrap(phaseTrue,[],2)*180/pi;
  phaseEst = unwrap(phaseEst*pi/180,[],2)*180/pi;
  k = round((phaseEst(:,2)-phaseTrue(:,2))/360);
  phaseEst = phaseEst - k(:,ones(1,size(phaseEst,2)))*360;

  % WH: compute and print RMSE between true and estimated freq. resp.
  magRMSE = 20*rms(log10(magTrue)-log10(magEst),2);
  magRMSE_all = 20*rms(log10(magTrue)-log10(magEst),"all");
  phaseRMSE = rms(phaseEst-phaseTrue,2);
  phaseRMSE_all = rms(phaseEst-phaseTrue,"all");
  fprintf(' - RMSE mag [dB] : %.2f (all) / ', magRMSE_all)
  for k = 1:length(magRMSE)
    fprintf('%.2f ', magRMSE(k));
  end
  fprintf('\n');
  fprintf(' - RMSE phs [deg]: %.2f (all) / ', phaseRMSE_all)
  for k = 1:length(magRMSE)
    fprintf('%.2f ', phaseRMSE(k));
  end
  fprintf('\n');

  % figure; clf;
  % subplot(1,2,1);
  % semilogx(W,20*log10(magTrue(:,:)),'linewidth',1.5); hold on
  % set(gca,'colororderindex',1);
  % semilogx(W,20*log10(abs(magEst(:,:))),'x');
  % semilogx([pi/Ts pi/Ts],ylim,'k:');
  % xlim([min(W) 1.1*max(W)]); grid on
  % xlabel('Analog frequency (rad/s)'); ylabel('Magnitude (dB)');
  % title(sprintf('Magnitude plots: Line = truth; x = est.; SOC = %d%%',100*SOC));
  %
  % subplot(1,2,2);
  % semilogx(W,phaseTrue,'linewidth',1.5); hold on
  % set(gca,'colororderindex',1);
  % semilogx(W,phaseEst,'x');
  % semilogx([pi/Ts pi/Ts],ylim,'k--');
  % xlim([min(W) 1.1*max(W)]); grid on
  % xlabel('Analog frequency (rad/s)'); ylabel('Phase (deg)');
  % title(sprintf('Phase plots: Line = truth; x = est.; SOC = %d%%',100*SOC));
  %
  % % WH: use thesisFormat if available.
  % if exist('thesisFormat', 'file')
  %   thesisFormat;
  % end

  close all % optional, but keeps total number of open plots small
  for mm = 1:size(magTrue,1)
    % if ~any(mm == [12, 18, 19, 20, 21])
    %     continue;
    % end

    figure; clf;
    subplot(1,2,1);
    semilogx(W,20*log10(magTrue(mm,:)),'linewidth',1.5); hold on
    set(gca,'colororderindex',1);
    semilogx(W,20*log10(abs(magEst(mm,:))),'x');
    semilogx([pi/Ts pi/Ts],ylim,'k:');
    xlim([min(W) 1.1*max(W)]); grid on
    xlabel('Analog frequency (rad/s)'); ylabel('Magnitude (dB)');
    title(sprintf('Magnitude plots: Line = truth; x = est.; SOC = %d%%',100*SOC));

    subplot(1,2,2);
    semilogx(W,phaseTrue(mm,:),'linewidth',1.5); hold on
    set(gca,'colororderindex',1);
    semilogx(W,phaseEst(mm,:),'x');
    semilogx([pi/Ts pi/Ts],ylim,'k--');
    xlim([min(W) 1.1*max(W)]); grid on
    xlabel('Analog frequency (rad/s)'); ylabel('Phase (deg)');
    title(sprintf('Phase plots: Line = truth; x = est.; SOC = %d%%',100*SOC));
    sgtitle(sprintf('TF # = %d, TF name = %s, loc = %d, kmax = %d',...
      mm,tfData.names{mm},tfData.xLoc(mm),kmaxlist(mm)));

    % WH: use thesisFormat if available.
    if exist('thesisFormat', 'file')
      try
        thesisFormat('PlotBoxPaddingInches', [0 0 0 0.3]);
      catch
        thesisFormat([0 0 0 0.3]);
      end % try
    end % if
  end % for

  drawnow;
end % if debug
end % xraHRA()

% -------------------------------------------------------------------------
% WH: Automatic determination of Kmax.
function [kmax, wc, tflistFlat] = autoKmax(tflist, wnyquist, cellData, Gc)
%AUTOKMAX Compute Kmax for each TF by finding the frequency at which
% the gain falls below the specified Gc (dB).

% First, split out TFs by x-location (one TF entry per x-location).
tflistFlat = {};
cursor = 1;
for ntf = 1:length(tflist)
  TFinfo = strsplit(tflist{ntf},{'tf','(s,',',%s)'});
  tfname = sprintf('tf%s',TFinfo{2});
  xlocs  = str2num(TFinfo{end-1})'; %#ok<ST2NM> (DO NOT USE STR2DOUBLE)
  for x = xlocs(:).'
    tflistFlat{cursor} = sprintf('%s(s,%d,%%s)',tfname,x); %#ok<AGROW>
    cursor = cursor + 1;
  end % for
end % for
numOut = length(tflistFlat);

% Identify frequency brackets in which each TF magnitude crosses Gc.
kmaxInit = 1;
wmax = wnyquist + 2*wnyquist*kmaxInit;
whigh = wmax*ones(numOut,1);
mask = true(numOut,1);
for i = 1:10
  freqResp = evalNormalizedTF(tflist,1j*wmax,cellData);
  magDB = 20*log10(abs(freqResp(:)));
  atCutoff = magDB <= Gc;
  whigh(mask) = wmax;
  if all(atCutoff)
    % All TFs below Gc at the maximum frequency.
    break;
  end % if
  mask = mask & ~atCutoff;

  % Increase wmax by one decade and re-evalulate frequency response.
  wmax = wnyquist + 2*wnyquist*kmaxInit*10^i;
end % for
wlow = whigh/10; % one decade below upper frequency

% Do binary search to find Gc frequencies (denoted wc).
wc = zeros(numOut,1);
for k = 1:numOut
  wl = wlow(k);
  wh = whigh(k);
  for i = 1:5
    wm = 10^mean(log10([wl wh]));
    freqResp = evalNormalizedTF(tflistFlat(k),1j*wm,cellData);
    gDB = 20*log10(abs(freqResp));
    if gDB <= Gc
      wh = wm;
    else
      wl = wm;
    end % if
  end % for i
  wc(k) = wm;
end % for numOut

% Compute kmax.
kmax = ceil((wc/wnyquist - 1)/2);
end % autoKmax()

% -------------------------------------------------------------------------
% WH: get TFs without D or dc gain.
function freqResp = evalNormalizedTF(tflist,s,cellData)
%EVALNORMALIZEDTF

if ~any(s==0)
  % Add zero frequency if not present.
  ss = [0 s(:).'];
end

[freqResp,hfGain] = evalTF(tflist,ss,cellData);
numOut = size(freqResp,1);

% Remove D term.
freqResp = (freqResp - diag(hfGain)*ones(numOut,length(ss))).';

% Set dc gain to 1.
dcGain = real(freqResp(ss==0,:)); % This is dcgain without D term
dcGain(dcGain==0) = 1;
dcGain(isinf(dcGain)) = 1;
dcGain(isnan(dcGain)) = 1;
freqResp(ss==0,:) = dcGain;
K = diag(dcGain);
freqResp = freqResp/K;

if ~any(s==0)
  % Remove zero frequency if not in requested frequency list.
  freqResp = freqResp(end,:);
end
end % evalNormalizedTF()

% -------------------------------------------------------------------------
% WH: C matrix optimization with fmincon: cost function and optimizer
% options.
function [fn, opts] = getCmatrixCostFunction(A,B,G,W,Ts)
% Cost function.
  function [J, H] = costC(c)
    J = 0;
    H = zeros(size(G));
    for kf = 1:length(W)
      Gk = G(kf);
      Mk = diag(1./(exp(1j*W(kf)*Ts)-diag(A)))*B;
      Hk = c*Mk;
      den = 1.0; %abs(Gk);
      J = J + (abs(Gk - Hk)./den).^2;
      H(kf) = Hk;
    end % for
  end % costC()
fn = @costC;

% Optimizer options.
opts = optimoptions( ...
  'fmincon', ...
  'Display','none', ...
  'FunctionTolerance',1e-12, ...
  'MaxFunctionEvaluations',100000 ...
  );
end % getCmatrixCostFunction()
