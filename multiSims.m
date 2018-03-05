% Object-oriented approach to simulating the control of the longitudinal
% dynamics of a very-flexible aircraft (VFA).
% 
% This script generates all the data required to run the simulations
% "VFA_ActOrder1.slx" and "VFA_ActOrder2.slx", which are
% relative-degree 2 and 3 MIMO VFA simulations, respectively.
% 
% The simulation uses a nonlinear VFA model as the plant, with controllers
% based on linearizations of the (nominal version of the) model at
% different dihedral angles. The adaptive controller is designed according
% to the method described in the paper "Adaptive Output-Feedback Control 
% for a Class of Multi-Input Multi-Output Plants with Arbitrary Relative 
% Degree" by Zheng (Max) Qu et al. These simulation files originate with
% Max's work, but have been continued by Ben Thomsen and other lab members
% since Max's graduation.
% 
% The main steps of the simulation process are:
%   1. Trim the nominal nonlinear VFA model (without actuator dynamics) at 
%       different dihedral angles (find the inputs and states such that 
%       states are not changing) and compute linearized state-space 
%       representation
%   2. Introduce uncertainties and actuator dynamics
%   3. For each trim point:
%       a. Augment linearized system with actuator dynamics and integral
%           error for command tracking
%       b. Add ficticious inputs to "square-up" the system - to provide
%           SPR properties necessary for adaptive control design
%       c. Compute baseline LQR feedback matrix (K), CRM feedback matrix 
%           (L), along with ouptut-mixing matrix (S) and transformed 
%           coordinates according to method in paper
%   4. Set initial conditions, choose adaptation rates and related gains
%       for simulation, define command trajectory, and run simulation
% 
% With first-order actuator model ("Relative-degree 2"): 
% States (x): Airspeed, Angle of Attack, Pitch Angle, Pitch Rate, Dihedral,
%             Dihedral Rate, Outer Aileron Deflection, Center Elevator 
%             Deflection, Dihedral Integral Error, Vert Accel Integral Error
%
% With second-order actuator model ("Relative-degree 3"):
% States (x): Airspeed, Angle of Attack, Pitch Angle, Pitch Rate, Dihedral, 
%             Dihedral Rate, Outer Aileron Deflection, Center Elevator 
%             Deflection, Outer Aileron Rate, Center Elevator Rate, 
%             Dihedral Integral Error, Vert Accel Integral Error
% 
% Outputs (y):      Pitch Rate, Dihedral Integral Error, 
%                   Vertical Accel Integral Error
%
% Goal is for outputs z to track command r (some notation uses z_cmd)
% Reg. Outputs (z): Dihedral Angle, Vertical Accel.
% 
% Inputs (u):       Outer Aileron, Center Elevator
% 
% The use of this simulation requires the Control System Toolbox and also
% Simulink Control Design
% 
% Relevant references are:
%       "Adaptive Control for a Class of Multi-Input Multi-Output Plants 
%           with Arbitrary Relative Degree"
%       
%       "Adaptive Output-Feedback Control and Applications to Very Flexible
%           Aircraft" (http://hdl.handle.net/1721.1/104223)
%
%       "Modeling for Control of Very Flexible Aircraft"
%           (https://doi.org/10.2514/6.2011-6202)
% 
%       "Squaring-Up Method in the Presence of Transmission Zeros"
%           (https://arxiv.org/abs/1310.1439v2)
% 
%       "Robust and Adaptive Control - with Aerospace Applications"
%           (https://doi.org/10.1007/978-1-4471-4396-3)

clear; clc;

opt.dataPath   = ['data', filesep]; % where to look for/store .mat files 
                                    % (relative to working directory)
opt.adaFlag    = true;              % use adaptive control?
opt.uncertFlag = true;              % uncertainty in dynamics?
opt.reComp     = true;              % trim, linearize, recompute controller
opt.pActOrder  = 1;                 % order of plant actuator dynamics (1 or 2)
opt.mActOrder  = 1;                 % order of modeled actuator dynamics
                                    % (1 or 2) and <= pActOrder
                                    
vfa1 = SimVFA(opt);  % initialize and setup the simulation
vfa1.runSim();       % run the simulation

%% Second Sim
opt = [];

opt.dataPath   = ['data', filesep]; % where to look for/store .mat files 
                                    % (relative to working directory)
opt.adaFlag    = true;              % use adaptive control?
opt.uncertFlag = true;              % uncertainty in dynamics?
opt.reComp     = true;              % trim, linearize, recompute controller
opt.pActOrder  = 2;                 % order of plant actuator dynamics (1 or 2)
opt.mActOrder  = 1;                 % order of modeled actuator dynamics
                                    % (1 or 2) and <= pActOrder

init_cond.tstart = 600;
ind1 = find(vfa1.simOutObj.t_sim>=init_cond.tstart,1);

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

%% Third Sim
opt = [];

opt.dataPath   = ['data', filesep]; % where to look for/store .mat files 
                                    % (relative to working directory)
opt.adaFlag    = true;              % use adaptive control?
opt.uncertFlag = true;              % uncertainty in dynamics?
opt.reComp     = true;              % trim, linearize, recompute controller
opt.pActOrder  = 2;                 % order of plant actuator dynamics (1 or 2)
opt.mActOrder  = 2;                 % order of modeled actuator dynamics
                                    % (1 or 2) and <= pActOrder

init_cond = [];
init_cond.tstart = 800;
ind2 = find(vfa2.simOutObj.t_sim>=init_cond.tstart,1);

init_cond.state_carry = zeros(7,1);
init_cond.state_carry(vfa2.simOpt.i_state_sel) = vfa2.simOutObj.xp(ind2,:);
init_cond.state_carry(3) = vfa2.trimPts.hinitial;
init_cond.xm_carry = vfa2.simOutObj.xm(1:length(vfa2.simOpt.i_state_sel), ind2);
init_cond.cmd_carry = [vfa2.simOutObj.r_cmd(:,ind2); vfa2.simOutObj.r_cmd_dot(:,ind2)];
init_cond.act_carry = vfa2.simOutObj.xact(:,ind2);

vfa3 = SimVFA(opt, init_cond);  % initialize and setup the simulation
vfa3.runSim();       % run the simulation

%% Data concatenation

tsim = vfa1.simOpt.tsim;
SOO1 = vfa1.simOutObj;
SOO2 = vfa2.simOutObj;
SOO3 = vfa3.simOutObj;

dihedral = (180/pi) * [SOO1.z(1,1:ind1), SOO2.z(1,1:ind2), SOO3.z(1,:)] + vfa1.simOpt.eta_nom;
vert_acc = [SOO1.z(2,1:ind1), SOO2.z(2,1:ind2), SOO3.z(2,:)];
out_ail  = (180/pi) * [SOO1.u_p(1,1:ind1), SOO2.u_p(1,1:ind2), SOO3.u_p(1,:)];
in_elev  = (180/pi) * [SOO1.u_p(2,1:ind1), SOO2.u_p(2,1:ind2), SOO3.u_p(2,:)];

u_ail  = (180/pi) * [SOO1.u_ad(1,1:ind1), SOO2.u_ad(1,1:ind2), SOO3.u_ad(1,:)];
u_elev = (180/pi) * [SOO1.u_ad(2,1:ind1), SOO2.u_ad(2,1:ind2), SOO3.u_ad(2,:)];

time = [SOO1.t_sim(1:ind1); SOO2.t_sim(1:ind2); SOO3.t_sim];
err_norm = [vecnorm(squeeze(SOO1.y(:,:,1:ind1)-SOO1.ym(:,:,1:ind1))), vecnorm(squeeze(SOO2.y(:,:,1:ind2)-SOO2.ym(:,:,1:ind2))), vecnorm(squeeze(SOO3.y-SOO3.ym))];
err_norm2 = [vecnorm(squeeze(SOO1.z(:,:,1:ind1)-SOO1.r_cmd(:,:,1:ind1))), vecnorm(squeeze(SOO2.z(:,:,1:ind2)-SOO2.r_cmd(:,:,1:ind2))), vecnorm(squeeze(SOO3.z-SOO3.r_cmd))];

%% Plotting: state and control

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

c2 = [0.466, 0.674, 0.188];

f1a = figure('Position',[1,1, 800,400]);
subplot(2,2,1)
plot(SOO1.t_sim, SOO1.r_cmd(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1.5)
hold on; grid on; 
plot(time, dihedral, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', c2)
xlim([0 tsim])
ylim([9 13])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Dihedral (deg)')
h=legend('Command', 'Output');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(2,2,2)
plot(SOO1.t_sim, SOO1.r_cmd(2,:), 'LineWidth', 1.5)
hold on; grid on; 
plot(time, vert_acc, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', c2)
xlim([0 tsim])
ylim([-3 3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Vertical Accel (ft/s^2)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(2,2,3)
plot(time, out_ail, 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]); grid on; hold on;
xlim([0 tsim])
ylim([0 3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Outer Aileron Command (deg)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(2,2,4)
plot(time, in_elev, 'LineWidth', 1.5, 'Color', [0, 0.447, 0.741]); grid on; hold on;
xlim([0 tsim])
ylim([-3 0])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Center Elevator Command (deg)')
xlabel('Time (s)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
% tightfig(f1a);


%% Plotting: error

f2 = figure('Position',[1,1, 800, 240]);
yyaxis left
plot(time, err_norm, 'LineWidth', 1.5); grid on; hold on;
xlim([0 tsim])
ylim([0 4e-3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
yyaxis right
plot(time, err_norm2, 'LineWidth', 1.5, 'Color', c2); grid on; hold on;
ylim([0 15])
title('Output Error Signals: $\|y(t)-y_m(t)\|_2$ and $\|z(t)-z_{cmd}(t)\|_2$', 'interpreter','latex')
xlabel('Time (s)')
h=legend('$\|y(t)-y_m(t)\|_2$', '$\|z(t)-z_{cmd}(t)\|_2$');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
% tightfig(f2);

%% Data concatenation for adaptive parameters

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

%% Plotting: adaptive parameters

figure('Position',[100,100, 800, 400]);
plot(time(1:ind1+ind2), norms_1, 'LineWidth', 1); grid on; hold on;
plot(time(ind1+ind2+1:end), norms_2, 'LineWidth', 1);
xlim([0 tsim]);
ylim([0 4]);
title('Normalized Learned Parameters', 'interpreter', 'latex')
xlabel('Time (s)')
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
h=legend('$\|\underline{\it{\Lambda}}\|$', '$\|\underline{\Psi}_1\|$', '$\|\underline{\Psi}_2\|$', '$\|\psi_{2}^1\|$', ...
         '$\|\underline{\it{\Lambda}}\|$', '$\|\underline{\Psi}_1\|$', '$\|\underline{\Psi}_2\|$', '$\|\psi_3^1\|$', '$\|\psi_3^2\|$', '$\|\underline{\Psi}_3\|$');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)


%% Bonus: VFA1 and VFA2 until loss of stability

SOO.t_sim(1:ind1) = vfa1.simOutObj.t_sim(1:ind1);
SOO.t_sim(ind1+1:length(vfa2.simOutObj.t_sim)+ind1) = vfa2.simOutObj.t_sim;

SOO.z(:,:,1:ind1) = vfa1.simOutObj.z(:,:,1:ind1);
SOO.z(:,:,ind1+1:ind1+length(vfa2.simOutObj.t_sim)) = vfa2.simOutObj.z(:,:,:);
SOO.u_p(:,:,1:ind1) = vfa1.simOutObj.u_p(:,:,1:ind1);
SOO.u_p(:,:,ind1+1:ind1+length(vfa2.simOutObj.t_sim)) = vfa2.simOutObj.u_p(:,:,:);

c2 = [0.466, 0.674, 0.188];

figure('Position',[1,1, 800, 400]);
subplot(2,2,1)
plot(SOO1.t_sim, SOO1.r_cmd(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1.5)
hold on; grid on; plot(SOO.t_sim, SOO.z(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1.5, 'Color', c2)
xlim([0 tsim])
ylim([9 13])
line([SOO.t_sim(ind1) SOO.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Dihedral (deg)')
h=legend('Command', 'Output');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(2,2,2)
plot(SOO1.t_sim, SOO1.r_cmd(2,:), 'LineWidth', 1.5)
hold on; grid on; plot(SOO.t_sim, SOO.z(2,:), 'LineWidth', 1.5, 'Color', c2)
xlim([0 tsim])
ylim([-3 3])
line([SOO.t_sim(ind1) SOO.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Vertical Accel (ft/s^2)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(2,2,3)
plot(SOO.t_sim, SOO.u_p(1,:)*180/pi, 'LineWidth', 1.5); grid on;
xlim([0 tsim])
ylim([0 3])
line([SOO.t_sim(ind1) SOO.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Outer Aileron (deg)')
xlabel('Time (s)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(2,2,4)
plot(SOO.t_sim, SOO.u_p(2,:)*180/pi, 'LineWidth', 1.5); grid on;
xlim([0 tsim])
ylim([-3 0])
line([SOO.t_sim(ind1) SOO.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Center Elevator (deg)')
xlabel('Time (s)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)


norm_lambda_ada = zeros(ind1 + length(vfa2.simOutObj.t_sim), 1);
norm_psi1_ada   = zeros(ind1 + length(vfa2.simOutObj.t_sim), 1);
norm_psi2_ada   = zeros(ind1 + length(vfa2.simOutObj.t_sim), 1); 
norm_psi21_ada = zeros(ind1 + length(vfa2.simOutObj.t_sim), 1);

for i=1:ind1
    norm_lambda_ada(i) = norm(SOO1.lambda_ada(:,:,i));
    norm_psi1_ada(i) = norm(SOO1.psi1_ada(:,:,i));
    norm_psi2_ada(i) = norm(SOO1.psi2_ada(:,:,i));
    norm_psi21_ada(i) = norm(SOO1.psi21_ada(:,:,i));
end

for i=1:length(vfa2.simOutObj.t_sim)
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
    
figure('Position',[100,100, 800, 240]);
plot(time(1:ind1+length(vfa2.simOutObj.t_sim)), norms_1, 'LineWidth', 1.5); grid on; hold on;
xlim([0 tsim]);
ylim([0 5]);
line([SOO.t_sim(ind1) SOO.t_sim(ind1)], ylim, 'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1);
title('Normalized Learned Parameters', 'interpreter', 'latex')
xlabel('Time (s)')
h=legend('$\|\underline{\it{\Lambda}}\|$', '$\|\underline{\Psi}_1\|$', '$\|\underline{\Psi}_2\|$', '$\|\psi_{2}^1\|$', ...
         '$\|\underline{\it{\Lambda}}\|$', '$\|\underline{\Psi}_1\|$', '$\|\underline{\Psi}_2\|$', '$\|\psi_3^1\|$', '$\|\psi_3^2\|$', '$\|\underline{\Psi}_3\|$');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)


err_norm_mis = [vecnorm(squeeze(SOO1.y(:,:,1:ind1)-SOO1.ym(:,:,1:ind1))), vecnorm(squeeze(SOO2.y(:,:,:)-SOO2.ym(:,:,:)))];
err_norm_mis2 = [vecnorm(squeeze(SOO1.z(:,:,1:ind1)-SOO1.r_cmd(:,:,1:ind1))), vecnorm(squeeze(SOO2.z(:,:,:)-SOO2.r_cmd(:,:,:)))];

f2 = figure('Position',[1,1, 800, 240]);
yyaxis left
plot(SOO.t_sim, err_norm_mis, 'LineWidth', 1.5); grid on; hold on;
xlim([0 tsim])
ylim([0 4e-3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off');
yyaxis right
plot(SOO.t_sim, err_norm_mis2, 'LineWidth', 1.5, 'Color', c2); grid on; hold on;
title('Output Error Signals: $\|y(t)-y_m(t)\|_2$ and $\|z(t)-z_{cmd}(t)\|_2$', 'interpreter','latex')
xlabel('Time (s)')
h=legend('$\|y(t)-y_m(t)\|_2$', '$\|z(t)-z_{cmd}(t)\|_2$');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)
% tightfig(f2);
