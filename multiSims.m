% Run the simulations to generate plots corresponding to: 
%   "Anomaly Response 3", "Anomaly Response 1", "Nom-3" in 
%   [Thomsen et al., IFAC CPHS 2018].
%
%   Anomaly Response w/o MRAC (RSLQR only) can be simulated by turning off
%   opt.adaFlag in VFA1, VFA2, VFA3 sections and commenting out adaptive
%   parameter carryover in VFA2 section (lines 44-48 as of 3/26/2018)

clear; clc;
opt.dataPath   = ['data', filesep]; % where to look for/store .mat files 
                                    % (relative to working directory)
                                    
%% VFA1
% Nominal operation (first-order actuators, first-order model) with
% uncertainty. RSLQR + MRAC. Start at t=0s.

opt.adaFlag    = true;              % use adaptive control?
opt.uncertFlag = true;              % uncertainty in dynamics?
opt.reComp     = true;              % trim, linearize, recompute controller
opt.pActOrder  = 1;                 % order of plant actuator dynamics (1 or 2)
opt.mActOrder  = 1;                 % order of modeled actuator dynamics
                                    % (1 or 2) and <= pActOrder
                                    
vfa1 = SimVFA(opt);  % initialize and setup the simulation
vfa1.runSim();       % run the simulation

%% VFA2: Second-order Actuators
% Anomalous operation (second-order actuators, first-order model) with
% uncertainty. RSLQR + MRAC. Start at t=600s.

opt.pActOrder  = 2; % change the order of the plant

init_cond.tstart = 600;
ind1 = find(vfa1.simOutObj.t_sim>=init_cond.tstart,1);

% carry over simulation state from VFA1 at t=tstart for continuity
init_cond.state_carry = zeros(7,1);
init_cond.state_carry(vfa1.simOpt.i_state_sel) = vfa1.simOutObj.xp(ind1,:);
init_cond.state_carry(3) = vfa1.trimPts.hinitial;
init_cond.xm_carry = vfa1.simOutObj.xm(1:length(vfa1.simOpt.i_state_sel), ind1);
init_cond.cmd_carry = [vfa1.simOutObj.r_cmd(:,ind1); vfa1.simOutObj.r_cmd_dot(:,ind1)];
init_cond.act_carry = vfa1.simOutObj.xact(:,ind1);

ada.lambda = vfa1.simOutObj.lambda_ada(:,:,ind1);
ada.psi1 = vfa1.simOutObj.psi1_ada(:,:,ind1);
ada.psi2 = vfa1.simOutObj.psi2_ada(:,:,ind1);
ada.psi21 = vfa1.simOutObj.psi21_ada(:,:,ind1);
init_cond.ada_carry = ada;

vfa2 = SimVFA(opt, init_cond);  % initialize and setup the simulation
vfa2.runSim();       % run the simulation

%% VFA3
% Anomalous operation (second-order actuators, second-order model) with
% uncertainty. RSLQR + MRAC. Start at t=800s.

opt.mActOrder  = 2;  % change control model of actuators to 2nd-order

init_cond = [];
init_cond.tstart = 800;
ind2 = find(vfa2.simOutObj.t_sim>=init_cond.tstart,1);

% carry over simulation state from VFA2 at t=tstart for continuity
init_cond.state_carry = zeros(7,1);
init_cond.state_carry(vfa2.simOpt.i_state_sel) = vfa2.simOutObj.xp(ind2,:);
init_cond.state_carry(3) = vfa2.trimPts.hinitial;
init_cond.xm_carry = vfa2.simOutObj.xm(1:length(vfa2.simOpt.i_state_sel), ind2);
init_cond.cmd_carry = [vfa2.simOutObj.r_cmd(:,ind2); vfa2.simOutObj.r_cmd_dot(:,ind2)];
init_cond.act_carry = vfa2.simOutObj.xact(:,ind2);

vfa3 = SimVFA(opt, init_cond);  % initialize and setup the simulation
vfa3.runSim();       % run the simulation

%% Concatenate signals

SOO1 = vfa1.simOutObj;
SOO2 = vfa2.simOutObj;
SOO3 = vfa3.simOutObj;

dihedral = (180/pi) * [SOO1.z(1,1:ind1), SOO2.z(1,1:ind2), SOO3.z(1,:)] + vfa1.simOpt.eta_nom;
vert_acc = [SOO1.z(2,1:ind1), SOO2.z(2,1:ind2), SOO3.z(2,:)];
out_ail  = (180/pi) * [SOO1.u_p(1,1:ind1), SOO2.u_p(1,1:ind2), SOO3.u_p(1,:)];
in_elev  = (180/pi) * [SOO1.u_p(2,1:ind1), SOO2.u_p(2,1:ind2), SOO3.u_p(2,:)];

% u_ail  = (180/pi) * [SOO1.u_ad(1,1:ind1), SOO2.u_ad(1,1:ind2), SOO3.u_ad(1,:)];
% u_elev = (180/pi) * [SOO1.u_ad(2,1:ind1), SOO2.u_ad(2,1:ind2), SOO3.u_ad(2,:)];

time = [SOO1.t_sim(1:ind1); SOO2.t_sim(1:ind2); SOO3.t_sim];
% model following error:
err_norm = [vecnorm(squeeze(SOO1.y(:,:,1:ind1)-SOO1.ym(:,:,1:ind1))), vecnorm(squeeze(SOO2.y(:,:,1:ind2)-SOO2.ym(:,:,1:ind2))), vecnorm(squeeze(SOO3.y-SOO3.ym))];
% tracking error:
err_norm2 = [vecnorm(squeeze(SOO1.z(:,:,1:ind1)-SOO1.r_cmd(:,:,1:ind1))), vecnorm(squeeze(SOO2.z(:,:,1:ind2)-SOO2.r_cmd(:,:,1:ind2))), vecnorm(squeeze(SOO3.z-SOO3.r_cmd))];


norm_lambda_ada = zeros(ind1 + ind2, 1);
norm_psi1_ada   = zeros(ind1 + ind2, 1);
norm_psi2_ada   = zeros(ind1 + ind2, 1); 
norm_psi21_ada = zeros(ind1 + ind2, 1);

for i=1:ind1
    norm_lambda_ada(i) = norm(SOO1.lambda_ada(:,:,i));
    norm_psi1_ada(i) = norm(SOO1.psi1_ada(:,:,i));
    norm_psi2_ada(i) = norm(SOO1.psi2_ada(:,:,i));
    norm_psi21_ada(i) = norm(SOO1.psi21_ada(:,:,i));
end

for i=1:ind2
    norm_lambda_ada(ind1+i) = norm(SOO2.lambda_ada(:,:,i));
    norm_psi1_ada(ind1+i) = norm(SOO2.psi1_ada(:,:,i));
    norm_psi2_ada(ind1+i) = norm(SOO2.psi2_ada(:,:,i));
    norm_psi21_ada(ind1+i) = norm(SOO2.psi21_ada(:,:,i));
end

norm_lambda_ada = norm_lambda_ada/norm_lambda_ada(ind1);
norm_psi1_ada   = norm_psi1_ada/norm_psi1_ada(ind1);
norm_psi2_ada   = norm_psi2_ada/norm_psi2_ada(ind1);
norm_psi21_ada  = norm_psi21_ada/norm_psi21_ada(ind1);

norms_1 = [norm_lambda_ada, norm_psi1_ada, norm_psi2_ada, norm_psi21_ada];
    
steps_2 = length(SOO3.t_sim);

norm_lambda_ada = zeros(steps_2, 1);
norm_psi1_ada   = zeros(steps_2, 1);
norm_psi2_ada   = zeros(steps_2, 1); 
norm_psi31_ada  = zeros(steps_2, 1);
norm_psi32_ada  = zeros(steps_2, 1);
norm_psi3_ada   = zeros(steps_2, 1);

for i=1:steps_2
    norm_lambda_ada(i) = norm(SOO3.lambda_ada(:,:,i));
    norm_psi1_ada(i) = norm(SOO3.psi1_ada(:,:,i));
    norm_psi2_ada(i) = norm(SOO3.psi2_ada(:,:,i));
    norm_psi31_ada(i) = norm(SOO3.psi31_ada(:,:,i));
    norm_psi32_ada(i) = norm(SOO3.psi32_ada(:,:,i));
    norm_psi3_ada(i) = norm(SOO3.psi3_ada(:,:,i));
end

norm_lambda_ada = norm_lambda_ada/norm_lambda_ada(end);
norm_psi1_ada   = norm_psi1_ada/norm_psi1_ada(end);
norm_psi2_ada   = norm_psi2_ada/norm_psi2_ada(end);
norm_psi31_ada  = norm_psi31_ada/norm_psi31_ada(end);
norm_psi32_ada  = norm_psi32_ada/norm_psi32_ada(end);
norm_psi3_ada   = norm_psi3_ada/norm_psi3_ada(end);

norms_2 = [norm_lambda_ada, norm_psi1_ada, norm_psi2_ada, ...
        norm_psi31_ada, norm_psi32_ada, norm_psi3_ada];
    
% set colors for plots
c1 = [0, 0.4470, 0.7410];
c2 = [0.466, 0.674, 0.188];

%% AR3 simulation plots: state + control

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

f1 = figure('Position',[1,1, 800, 640]);
subplot(4,1,1)
plot(SOO1.t_sim, SOO1.r_cmd(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1.5)
hold on; grid on; 
plot(time, dihedral, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', c2)
xlim([0 1500])
ylim([9 13])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(end) SOO2.t_sim(end)],ylim,'Color',[0.6 0.6 0.6],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
title('Dihedral (deg)','Interpreter','Latex')
h=legend('Command', 'Output');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(4,1,2)
plot(SOO1.t_sim, SOO1.r_cmd(2,:), 'LineWidth', 1.5)
hold on; grid on; 
plot(time, vert_acc, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', c2)
xlim([0 1500])
ylim([-3 3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(end) SOO2.t_sim(end)],ylim,'Color',[0.6 0.6 0.6],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
title('Vertical Accel (ft/s$^2$)','Interpreter','Latex')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(4,1,3)
plot(time, out_ail, 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]); grid on; hold on;
xlim([0 1500])
ylim([0 3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(end) SOO2.t_sim(end)],ylim,'Color',[0.6 0.6 0.6],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
title('Outer Aileron (deg)','Interpreter','Latex')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(4,1,4)
plot(time, in_elev, 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]); grid on; hold on;
xlim([0 1500])
ylim([-3 0])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(end) SOO2.t_sim(end)],ylim,'Color',[0.6 0.6 0.6],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
title('Center Elevator (deg)','Interpreter','Latex')
xlabel('Time (s)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
tightfig(f1);

%% AR3 simulation plots: errors

f2 = figure('Position',[1,1, 800, 240]);
set(f2,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(time, err_norm, 'LineWidth', 1.5); grid on; hold on;
xlim([0 1500])
ylim([0 4e-3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(end) SOO2.t_sim(end)],ylim,'Color',[0.6 0.6 0.6],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
yyaxis right
plot(time, err_norm2, 'LineWidth', 1.5, 'Color', c2); grid on; hold on;
ylim([0 15])
title('Output Error Signals: $\|y(t)-y_m(t)\|_2$ and $\|z(t)-z_{cmd}(t)\|_2$', 'interpreter','latex')
xlabel('Time (s)')
h=legend('$\|y(t)-y_m(t)\|_2$', '$\|z(t)-z_{cmd}(t)\|_2$');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
tightfig(f2);

%% AR3 simulation plots: adaptive parameters

f3 = figure('Position',[100,100, 800, 240]);
plot(time(1:ind1+ind2), norms_1, 'LineWidth', 1); grid on; hold on;
plot(time(ind1+ind2+1:end), norms_2, 'LineWidth', 1);
xlim([0 1500]);
ylim([0 4]);
title('Normalized Learned Parameters', 'interpreter', 'latex')
xlabel('Time (s)')
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
line([SOO2.t_sim(end) SOO2.t_sim(end)],ylim,'Color',[0.6 0.6 0.6],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
tightfig(f3);

%% AR1 simulation plots
% "Passive" anomaly response, showing VFA1 and then VFA2 until loss of
% stability -- just hacked together by copying and pasting from AR3 stuff

AR1.t_sim(1:ind1) = vfa1.simOutObj.t_sim(1:ind1);
AR1.t_sim(ind1+1:length(vfa2.simOutObj.t_sim)+ind1) = vfa2.simOutObj.t_sim;

AR1.z(:,:,1:ind1) = vfa1.simOutObj.z(:,:,1:ind1);
AR1.z(:,:,ind1+1:ind1+length(vfa2.simOutObj.t_sim)) = vfa2.simOutObj.z(:,:,:);
AR1.u_p(:,:,1:ind1) = vfa1.simOutObj.u_p(:,:,1:ind1);
AR1.u_p(:,:,ind1+1:ind1+length(vfa2.simOutObj.t_sim)) = vfa2.simOutObj.u_p(:,:,:);

AR1.norm_lambda_ada = zeros(ind1 + length(vfa2.simOutObj.t_sim), 1);
AR1.norm_psi1_ada   = zeros(ind1 + length(vfa2.simOutObj.t_sim), 1);
AR1.norm_psi2_ada   = zeros(ind1 + length(vfa2.simOutObj.t_sim), 1); 
AR1.norm_psi21_ada = zeros(ind1 + length(vfa2.simOutObj.t_sim), 1);

for i=1:ind1
    AR1.norm_lambda_ada(i) = norm(SOO1.lambda_ada(:,:,i));
    AR1.norm_psi1_ada(i) = norm(SOO1.psi1_ada(:,:,i));
    AR1.norm_psi2_ada(i) = norm(SOO1.psi2_ada(:,:,i));
    AR1.norm_psi21_ada(i) = norm(SOO1.psi21_ada(:,:,i));
end

for i=1:length(vfa2.simOutObj.t_sim)
    AR1.norm_lambda_ada(ind1+i) = norm(SOO2.lambda_ada(:,:,i));
    AR1.norm_psi1_ada(ind1+i) = norm(SOO2.psi1_ada(:,:,i));
    AR1.norm_psi2_ada(ind1+i) = norm(SOO2.psi2_ada(:,:,i));
    AR1.norm_psi21_ada(ind1+i) = norm(SOO2.psi21_ada(:,:,i));
end

AR1.norm_lambda_ada = AR1.norm_lambda_ada/AR1.norm_lambda_ada(ind1);
AR1.norm_psi1_ada   = AR1.norm_psi1_ada/AR1.norm_psi1_ada(ind1);
AR1.norm_psi2_ada   = AR1.norm_psi2_ada/AR1.norm_psi2_ada(ind1);
AR1.norm_psi21_ada  = AR1.norm_psi21_ada/AR1.norm_psi21_ada(ind1);

AR1.norms_1 = [AR1.norm_lambda_ada, AR1.norm_psi1_ada, AR1.norm_psi2_ada, AR1.norm_psi21_ada];

% model-following error:
AR1.err_norm = [vecnorm(squeeze(SOO1.y(:,:,1:ind1)-SOO1.ym(:,:,1:ind1))), vecnorm(squeeze(SOO2.y(:,:,:)-SOO2.ym(:,:,:)))];
% tracking error:
AR1.err_norm2 = [vecnorm(squeeze(SOO1.z(:,:,1:ind1)-SOO1.r_cmd(:,:,1:ind1))), vecnorm(squeeze(SOO2.z(:,:,:)-SOO2.r_cmd(:,:,:)))];

% AR1 simulation plots: state + control:
f1_ar1 = figure('Position',[1,1, 800, 640]);
subplot(4,1,1)
plot(SOO1.t_sim, SOO1.r_cmd(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1.5)
hold on; grid on; plot(AR1.t_sim, AR1.z(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1.5, 'Color', c2)
xlim([0 1500])
ylim([9 13])
line([AR1.t_sim(ind1) AR1.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Dihedral (deg)','Interpreter','Latex')
h=legend('Command', 'Output');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(4,1,2)
plot(SOO1.t_sim, SOO1.r_cmd(2,:), 'LineWidth', 1.5)
hold on; grid on; plot(AR1.t_sim, AR1.z(2,:), 'LineWidth', 1.5, 'Color', c2)
xlim([0 1500])
ylim([-3 3])
line([AR1.t_sim(ind1) AR1.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Vertical Accel (ft/s$^2$)','Interpreter','Latex')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(4,1,3)
plot(AR1.t_sim, AR1.u_p(1,:)*180/pi, 'LineWidth', 1.5); grid on;
xlim([0 1500])
ylim([0 3])
line([AR1.t_sim(ind1) AR1.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Outer Aileron (deg)','Interpreter','Latex')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(4,1,4)
plot(AR1.t_sim, AR1.u_p(2,:)*180/pi, 'LineWidth', 1.5); grid on;
xlim([0 1500])
ylim([-3 0])
line([AR1.t_sim(ind1) AR1.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Center Elevator (deg)','Interpreter','Latex')
xlabel('Time (s)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
tightfig(f1_ar1);

% AR1 simulation plots: errors
f2_ar1 = figure('Position',[1,1, 800, 240]);
set(f2_ar1,'defaultAxesColorOrder',[c1; c2]);
yyaxis left
plot(AR1.t_sim, AR1.err_norm, 'LineWidth', 1.5); grid on; hold on;
xlim([0 1500])
ylim([0 4e-3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([vfa2.simOutObj.t_sim(end) vfa2.simOutObj.t_sim(end)], ylim, 'Color',[0.3 0.3 0.3],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
yyaxis right
ylim([0 15])
plot(AR1.t_sim, AR1.err_norm2, 'LineWidth', 1.5, 'Color', c2); grid on; hold on;
title('Output Error Signals: $\|y(t)-y_m(t)\|_2$ and $\|z(t)-z_{cmd}(t)\|_2$', 'interpreter','latex')
xlabel('Time (s)')
h=legend('$\|y(t)-y_m(t)\|_2$', '$\|z(t)-z_{cmd}(t)\|_2$');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
tightfig(f2_ar1);

% AR1 simulation plots: adaptive parameters
f3_ar1 = figure('Position',[100,100, 800, 240]);
plot(AR1.t_sim, AR1.norms_1, 'LineWidth', 1.5); grid on; hold on;
xlim([0 1500]);
ylim([0 5]);
line([AR1.t_sim(ind1) AR1.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([vfa2.simOutObj.t_sim(end) vfa2.simOutObj.t_sim(end)], ylim, 'Color',[0.3 0.3 0.3],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
title('Normalized Learned Parameters', 'interpreter', 'latex')
xlabel('Time (s)')
h=legend('$\|\underline{\it{\Lambda}}\|$', '$\|\underline{\Psi}_1\|$', '$\|\underline{\Psi}_2\|$', '$\|\psi_{2}^1\|$');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
tightfig(f3_ar1);


%% Nom-3 simulation plots
% Nominal operation with RSLQR+MRAC, no anomalies.

% Nom-3 simulation plots: state + control:
f1_nom3 = figure('Position',[1,1, 800, 640]);
subplot(4,1,1)
plot(SOO1.t_sim, SOO1.r_cmd(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1.5)
hold on; grid on; plot(SOO1.t_sim, SOO1.z(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1.5, 'Color', c2)
xlim([0 1500])
ylim([9 13])
title('Dihedral (deg)','Interpreter','Latex')
h=legend('Command', 'Output');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(4,1,2)
plot(SOO1.t_sim, SOO1.r_cmd(2,:), 'LineWidth', 1.5)
hold on; grid on; plot(SOO1.t_sim, SOO1.z(2,:), 'LineWidth', 1.5, 'Color', c2)
xlim([0 1500])
ylim([-3 3])
title('Vertical Accel (ft/s$^2$)','Interpreter','Latex')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(4,1,3)
plot(SOO1.t_sim, SOO1.u_p(1,:)*180/pi, 'LineWidth', 1.5); grid on;
xlim([0 1500])
ylim([0 3])
title('Outer Aileron (deg)','Interpreter','Latex')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(4,1,4)
plot(SOO1.t_sim, SOO1.u_p(2,:)*180/pi, 'LineWidth', 1.5); grid on;
xlim([0 1500])
ylim([-3 0])
title('Center Elevator (deg)','Interpreter','Latex')
xlabel('Time (s)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
tightfig(f1_nom3);
