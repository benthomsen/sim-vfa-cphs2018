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
opt.tstart = 0;
                                    
vfa1 = SimVFA(opt);  % initialize and setup the simulation
vfa1.runSim();       % run the simulation
% vfa1.plotSim();      % plot output from the simulation

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

opt.tstart = 400;
ind1 = find(vfa1.simOutObj.t_sim>=opt.tstart,1);

opt.state_carry = zeros(7,1);
opt.state_carry(vfa1.simOpt.i_state_sel) = vfa1.simOutObj.xp(ind1,:);
opt.state_carry(3) = vfa1.trimPts.hinitial;
opt.xm_carry = vfa1.simOutObj.xm(1:length(vfa1.simOpt.i_state_sel), ind1);
opt.cmd_carry = [vfa1.simOutObj.r_cmd(:,ind1); vfa1.simOutObj.r_cmd_dot(:,ind1)];
opt.act_carry = vfa1.simOutObj.xact(:,ind1);

ada.lambda = vfa1.simOutObj.lambda_ada(:,:,ind1);
ada.psi1 = vfa1.simOutObj.psi1_ada(:,:,ind1);
ada.psi2 = vfa1.simOutObj.psi2_ada(:,:,ind1);
ada.psi21 = vfa1.simOutObj.psi21_ada(:,:,ind1);
opt.ada_carry = ada;

vfa2 = SimVFA(opt);  % initialize and setup the simulation
vfa2.runSim();       % run the simulation
% vfa2.plotSim();      % plot output from the simulation

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

opt.tstart = 550;
ind2 = find(vfa2.simOutObj.t_sim>=opt.tstart,1);

opt.state_carry = zeros(7,1);
opt.state_carry(vfa2.simOpt.i_state_sel) = vfa2.simOutObj.xp(ind2,:);
opt.state_carry(3) = vfa2.trimPts.hinitial;
opt.xm_carry = vfa2.simOutObj.xm(1:length(vfa2.simOpt.i_state_sel), ind2);
opt.cmd_carry = [vfa2.simOutObj.r_cmd(:,ind2); vfa2.simOutObj.r_cmd_dot(:,ind2)];
opt.act_carry = vfa2.simOutObj.xact(:,ind2);

vfa3 = SimVFA(opt);  % initialize and setup the simulation
vfa3.runSim();       % run the simulation
% vfa3.plotSim();      % plot output from the simulation

%% Plotting 

tsim = vfa1.simOpt.tsim;
SOO1 = vfa1.simOutObj;
SOO2 = vfa2.simOutObj;
SOO3 = vfa3.simOutObj;

figure('Position',[1,1, 1000, 400]);
subplot(2,2,1)
plot(SOO1.t_sim, SOO1.r_cmd(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1)
hold on; grid on; 
plot(SOO1.t_sim(1:ind1), SOO1.z(1,1:ind1)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1, 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098])
plot(SOO2.t_sim(1:ind2), SOO2.z(1,1:ind2)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1, 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098])
plot(SOO3.t_sim, SOO3.z(1,:)*180/pi + vfa1.simOpt.eta_nom, 'LineWidth', 1, 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098])
xlim([0 tsim])
ylim([9 13])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle',':', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle',':', 'LineWidth', 1);
title('Dihedral (deg)')
h=legend('Command', 'Output');
set(h,'fontsize',vfa1.pltOpt.legfontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname,'Interpreter','Latex','Location','SouthEast'); legend('boxoff')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(2,2,2)
plot(SOO1.t_sim, SOO1.r_cmd(2,:), 'LineWidth', 1)
hold on; grid on; 
plot(SOO1.t_sim(1:ind1), SOO1.z(2,1:ind1), 'LineWidth', 1, 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098])
plot(SOO2.t_sim(1:ind2), SOO2.z(2,1:ind2), 'LineWidth', 1, 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098])
plot(SOO3.t_sim, SOO3.z(2,:), 'LineWidth', 1, 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098])
xlim([0 tsim])
ylim([-1.5 1.5])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle',':', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle',':', 'LineWidth', 1);
title('Vertical Accel (ft/s^2)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(2,2,3)
plot(SOO1.t_sim(1:ind1), SOO1.u_p(1,1:ind1)*180/pi, 'LineWidth', 1); grid on; hold on;
plot(SOO2.t_sim(1:ind2), SOO2.u_p(1,1:ind2)*180/pi, 'LineWidth', 1);
plot(SOO3.t_sim, SOO3.u_p(1,:)*180/pi, 'LineWidth', 1);
xlim([0 tsim])
ylim([0 3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle',':', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle',':', 'LineWidth', 1);
title('Outer Aileron (deg)')
xlabel('Time (s)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)

subplot(2,2,4)
plot(SOO1.t_sim(1:ind1), SOO1.u_p(2,1:ind1)*180/pi, 'LineWidth', 1); grid on; hold on;
plot(SOO2.t_sim(1:ind2), SOO2.u_p(2,1:ind2)*180/pi, 'LineWidth', 1);
plot(SOO3.t_sim, SOO3.u_p(2,:)*180/pi, 'LineWidth', 1);
xlim([0 tsim])
ylim([0 3])
line([SOO1.t_sim(ind1) SOO1.t_sim(ind1)],ylim,'Color',[0 0 0],'LineStyle',':', 'LineWidth', 1);
line([SOO2.t_sim(ind2) SOO2.t_sim(ind2)],ylim,'Color',[0 0 0],'LineStyle',':', 'LineWidth', 1);
title('Center Elevator (deg)')
xlabel('Time (s)')
set(gca,'fontsize',vfa1.pltOpt.fontsize,'fontweight',vfa1.pltOpt.weight,'fontname',vfa1.pltOpt.fontname)