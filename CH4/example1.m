% This example reproduces the results in section 4.2 of the book. It uses
% the presented method to convert a continuous-time frequency response to a
% discrete-time frequency response.
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

clearvars; close all; clc;
if(~isdeployed),cd(fileparts(which(mfilename)));end
addpath(genpath('../'));

% first example from this book subsection:
% define continuous-time system, sample period, continuous-time frequency
% range, and compute the frequency response and "D" term
csys = tf([.1 1 1],[1 3 1]); % continuous-time transfer function
Ts = 0.1;                    % sampling period (10 Hz frequency)
w = [0,logspace(-6,6,1e4)];  % continuous-time frequency range
Gc(1,:) = freqresp(csys,w);  % continuous-time frquency response
Dss = ss(csys); D = Dss.D;   % compute infinite-frequency gain "D"

% compute the discrete-time frequency response using proposed general
% method and MATLAB's built-in polynomial method c2d
wnyquist = pi/Ts; % the Nyquist sampling frequency
W = [0,logspace(-5,0,50)]*wnyquist; % 51 frequencies b/w 0, Nyquist rate
kmax = 5; % number of aliases to consider +/- of "0" alias

Hd = c2dGeneral(w,Gc,kmax,W,Ts,D); % compute using proposed general method
Gd = squeeze(freqresp(c2d(csys,Ts),W)); % compute using polynomial method

% plot magnitude responses
figure(1); clf;
semilogx(w,20*log10(abs(Gc)),W,20*log10(abs(Gd))); hold on
set(gca,'colororderindex',2); 
semilogx(W,20*log10(abs(Hd)),'.','markersize',10);

xlim([min(W(2:end)),max(W)]); set(gca,'xtick',10.^(-3:1));
legend('Continuous-time system','Actual discrete-time system',...
  'Estimated discrete-time system','location','southwest');
xlabel('Frequency (rad s^{-1})'); ylabel('Magnitude (dB)');
title('Magnitude response of polynomial system'); grid on

% phase
figure(2); clf;
semilogx(w,angle(Gc)*180/pi,W,angle(Gd)*180/pi); hold on
set(gca,'colororderindex',2); 
semilogx(W,angle(Hd)*180/pi,'.','markersize',10);
xlim([min(W(2:end)),max(W)]); set(gca,'xtick',10.^(-3:1));
legend('Continuous-time system','Actual discrete-time system',...
  'Estimated discrete-time system','location','southwest');
xlabel('Frequency (rad s^{-1})'); ylabel('Phase (deg)');
title('Phase response of polynomial system'); grid on

% second example from this subsection
% define continuous-time system, sample period, continuous-time frequency
% range, and compute the frequency response and "D" term
csys = tf([1 10 1],[10 .1 0]); % continuous-time transfer function
Ts = 0.1;                    % sampling period (10 Hz frequency)
w = [0,logspace(-6,6,1e4)];  % continuous-time frequency range
Gc(1,:) = freqresp(csys,w);  % continuous-time frquency response
Dss = ss(csys); D = Dss.D;   % compute infinite-frequency gain "D"

% compute the discrete-time frequency response using proposed general
% methd and MATLAB's built-in polynomial method c2d
wnyquist = pi/Ts; % the Nyquist sampling frequency
W = [0,logspace(-5,0,50)]*wnyquist; % 51 frequencies b/w 0, Nyquist rate
kmax = 15; % number of aliases to consider +/- of "0" alias

Hd = c2dGeneral(w,Gc,kmax,W,Ts,D); % compute using proposed general method
warning('off','Control:analysis:InfiniteFreqResp');
Gd = squeeze(freqresp(c2d(csys,Ts),W)); % compute using polynomial method
warning('on','Control:analysis:InfiniteFreqResp');

% plot magnitude responses
figure(3); clf;
semilogx(w,20*log10(abs(Gc)),W,20*log10(abs(Gd))); hold on
set(gca,'colororderindex',2); 
semilogx(W,20*log10(abs(Hd)),'.','markersize',10);
xlim([min(W(2:end)),max(W)]); set(gca,'xtick',10.^(-3:1));
legend('Continuous-time system','Actual discrete-time system',...
  'Estimated discrete-time system','location','northeast');
xlabel('Frequency (rad s^{-1})'); ylabel('Magnitude (dB)');
title('Magnitude response of polynomial system'); grid on

% phase
figure(4); clf;
semilogx(w,angle(Gc)*180/pi,W,angle(Gd)*180/pi); hold on
set(gca,'colororderindex',2); 
semilogx(W,angle(Hd)*180/pi,'.','markersize',10);
xlim([min(W(2:end)),max(W)]); set(gca,'xtick',10.^(-3:1));
legend('Continuous-time system','Actual discrete-time system',...
  'Estimated discrete-time system','location','northwest');
xlabel('Frequency (rad s^{-1})'); ylabel('Phase (deg)');
title('Phase response of polynomial system'); grid on

% this is the function that does the actual conversion
function Hd = c2dGeneral(w,G,kmax,W,Ts,Dterm)
  % Compute continuous-time frequencies needed for method
  M = (ones(2*kmax+1,1)*W) + 2*pi/Ts*((-kmax:kmax)'*ones(1,length(W)));
  % Compute input frequency responses for negative frequencies 
  % Note that G(-w) = cong(G(w))... and don't repeat w=0 
  [winterp,ind1] = sort([-w(w~=0),w]); % negative frequencies
  Ginterp = [conj(G(:,w~=0)),G];       % frequency responses
  % Compute input frequency responses at all needed frequencies in "M"
  freqResp = interp1(winterp,Ginterp(:,ind1).',M(:)).';
  numOut = size(freqResp,1);
  % remove D term for now; add back into final model
  freqResp = (freqResp - diag(Dterm)*ones(numOut,length(M(:)))).';
  % Add together frequency aliases to create discrete-time freq resp
  Geval = [];
  for k = 1:numOut % glue together responses horizontally 
    newData = freqResp(:,k).*(1-exp(-1j*M(:)*Ts))./(1j*M(:));
    newData(M(:)==0) = freqResp(M(:)==0,k)*Ts;
    Geval = [Geval, reshape(newData.',size(M))]; %#ok<AGROW>
  end % Geval has rows for -kmax...kmax, repeated numOut times horizontally
  Gdest = (1/Ts)*sum(Geval);
  Gdest = reshape(Gdest.',length(W),numOut).'; % reshape outputs
  % Add D term back into true freq response for plotting purposes
  Hd = Gdest + diag(Dterm)*ones(numOut,length(W));
end  
