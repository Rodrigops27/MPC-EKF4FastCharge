function plot4Figures(time, u_store, voltage_store, SOC_store, phise_store, mpcData)
% plot4Figures.m
% 1: Current, 2: Voltage, 3: SOC  (same as plotFigures)
% 4: Side-reaction overpotential (Phise / eta)
%
% Inputs:
%   time          : N×1 time vector [s]
%   u_store       : N×1 current [A]
%   voltage_store : N×1 ROM voltage [V]
%   SOC_store     : 1×N or N×1 SOC trace (fraction 0–1 or %)
%   phise_store   : N×1 side-reaction overpotential [V] (your logged signal)
%   mpcData       : struct with .const fields for limits (u/v/phi), .ref (target SOC)

% --- make vectors consistent length ---
t   = time(:);
I   = u_store(:);
V   = voltage_store(:);
SOC = SOC_store(:);
Phi = phise_store(:);

n = min([numel(t), numel(I), numel(V), numel(SOC), numel(Phi)]);
t   = t(1:n);  I = I(1:n);  V = V(1:n);  SOC = SOC(1:n);  Phi = Phi(1:n);

% --- figure and axes ---
hFig = figure('Name','MPC Fast-Charge (4×)','Color','w','Position',[600 120 680 820]);
ax = gobjects(4,1);

% ---- 1) Current ----
ax(1) = subplot(4,1,1);
plot(t, I, 'LineWidth', 2, 'DisplayName','Charging current'); hold on;
if isfield(mpcData.const,'u_max'), yline(mpcData.const.u_max,'--r','LineWidth',1.5,'DisplayName','u_{max}'); end
if isfield(mpcData.const,'u_min'), yline(mpcData.const.u_min,'--r','LineWidth',1.5,'DisplayName','u_{min}'); end
grid on; title('Charging Current','FontSize',14);
ylabel('Current [A]','FontSize',12);
legend('Location','northeast');
xlim([t(1) t(end)]);

% ---- 2) Voltage ----
ax(2) = subplot(4,1,2);
plot(t, V, 'LineWidth', 2, 'DisplayName','Cell voltage'); hold on;
if isfield(mpcData.const,'v_max'), yline(mpcData.const.v_max,'--r','LineWidth',1.5,'DisplayName','v_{max}'); end
if isfield(mpcData.const,'v_min'), yline(mpcData.const.v_min,'--r','LineWidth',1.5,'DisplayName','v_{min}'); end
grid on; title('Voltage','FontSize',14);
ylabel('Voltage [V]','FontSize',12);
legend('Location','northeast');
xlim([t(1) t(end)]);

% ---- 3) SOC ----
ax(3) = subplot(4,1,3);
plot(t, SOC, 'LineWidth', 2, 'DisplayName','SOC'); hold on;
if isfield(mpcData,'ref') && ~isempty(mpcData.ref)
    yline(mpcData.ref/100,'--r','LineWidth',1.5,'DisplayName','Target');
end
grid on; title('SOC','FontSize',14);
xlabel(''); % bottom axes will carry the xlabel
if max(SOC) <= 1.5
    ylabel('SOC [–]','FontSize',12); ylim([0 1]);
else
    ylabel('SOC [%]','FontSize',12); ylim([0 100]);
end
legend('Location','northeast');
xlim([t(1) t(end)]);

% ---- 4) Side-reaction overpotential (Phise / eta) ----
ax(4) = subplot(4,1,4);
plot(t, Phi, 'LineWidth', 2, 'DisplayName','\phi_{se} / \eta'); hold on;
if isfield(mpcData.const,'phise_min')
    yline(mpcData.const.phise_min,'--r','LineWidth',1.5,'DisplayName','\phi_{min}');
end
yline(0,'k-','LineWidth',1); % reference line
grid on; title('Side-reaction Overpotential','FontSize',14);
ylabel('Overpot. [V]','FontSize',12);
xlabel('Time [s]','FontSize',12);
legend('Location','northeast');
xlim([t(1) t(end)]);

% ---- Link x-axes for synchronized zoom/pan ----
linkaxes(ax,'x');

end
