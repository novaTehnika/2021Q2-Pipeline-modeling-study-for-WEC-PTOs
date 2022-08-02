% run_PTOpipeline.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 03/29/2021
%
% PURPOSE/DESCRIPTION:
% This script 
%
% FILE DEPENDENCY:
% sim_PTOpipeline.m
% pLmodelParamSetup.m
% parameters_PTOpipelineStudy.m
%
% UPDATES:
% 03/29/2021 - Creation. 
% 03/31/2021 - Added tolerance parameters for ODE solvers to par structure.
% 03/30/2021 - Revised to use external function 
% ("parameters_PTOpipelineStudy.m") to specifiy parameters, which includes
% design case switching.
% 05/18/2021 - Updated model versions implimenting unsteady 
%           friction in clear cmd (2nd executable line).
%
% Copyright (C) 2022  Jeremy W. Simmons II
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clear pumpFlow pipeline pipelineMOC
% clc

%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Simulation Parameters
par.tend = 600; %[s] time span of simulation

par.odeSolverRelTol = 1e-6; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-6; % Abs. error tolerance parameter for ODE solver

rng(2); % set the seed for the random number generator
    
%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify simulation parameters
designCase = 5;
designSubCase = 1;
pLmodel = 10;
par = parameters_PTOpipeline(par,designCase,designSubCase,pLmodel);

%%%%%%%%%%%%
n_seg = 200;
par = pLmodelParamSetup(par,n_seg);
 %%%%%%%%%%%%%%%
 
% run simulation
tic
out = sim_PTOpipeline(par,1);
toc

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_mean = 1e-6*mean(out.p_hout)
PP_load = 1e-3*mean(out.q_p.*out.p_hout)
% q_mean = out.par.Xq/12
% mean(out.q_p)

vPL_mean = mean(out.q_hout/out.par.A_line)

if sum(par.pLmodel == [8 9 10 11])
    PP_HP = out.PPmean_HP_fric
    PP_HP_sFric = out.PPmean_HP_sfric
    PP_HP_uFric = out.PPmean_HP_ufric
end
figure
hold on
plot(out.t,out.p_hin)
plot(out.t,out.p_hout)
plot(out.t,out.p_lout)
legend('p_{hin}','p_{hout}','p_{lout}')
title(['Design Case = ',num2str(designCase),'.',num2str(designSubCase),', PL Model = ',num2str(pLmodel)])

%%
if 0
    
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = .5;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 3;

%%% at higher pressures
ax1 = subplot(n_plots,1,1);
ax1.FontName = 'times';
ax1.FontSize = fontSize-1;

% yyaxis left
hold on
grid on
plot(out.t-50,1e-6*out.p_hin,'k-','LineWidth',lineWidth)
plot(out.t-50,1e-6*out.p_hout,'k-.','LineWidth',lineWidth)
ylabel('Pressure (MPa)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

yyaxis right
hold on
grid on
plot(out.t-50,1e3*out.q_hin,'r-','LineWidth',lineWidth)
plot(out.t-50,1e3*out.q_hout,'r--','LineWidth',lineWidth)
ylabel('Flow rate (L/s)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

leg = legend('p_{h,in}','p_{h,out}','q_{h,in}','q_{h,out}');
leg.FontSize = fontSize-1;
leg.FontName = 'Cambria Math';

xlabel('Time (s)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
% title('High-Pressure Pipeline Pressure Dynamics',...
%     'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

%%% at lower pressures
ax2 = subplot(n_plots,1,2);
ax2.FontName = 'times';
ax2.FontSize = fontSize-1;

% yyaxis left
hold on
grid on
plot(out.t-50,1e-6*out.p_lin,'k-.','LineWidth',lineWidth)
plot(out.t-50,1e-6*out.p_lout,'k-','LineWidth',lineWidth)
ylabel('Pressure (MPa)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

yyaxis right
hold on
grid on
plot(out.t-50,1e3*out.q_lin,'r-','LineWidth',lineWidth)
plot(out.t-50,1e3*out.q_lout,'k--','LineWidth',lineWidth)
ylabel('Flow rate (L/s)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

leg = legend('p_{l,in}','p_{l,out}','q_{l,in}','q_{l,out}');
leg.FontSize = fontSize-1;
leg.FontName = 'Cambria Math';


xlabel('Time (s)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
% title('Low-Pressure Pipeline Pressure Dynamics',...
%     'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')


%%% at higher pressures
ax3 = subplot(n_plots,1,3);
ax3.FontName = 'times';
ax3.FontSize = fontSize-1;

% yyaxis left
hold on
grid on
plot(out.t-50,1e3*out.q_p,'k-','LineWidth',lineWidth)
plot(out.t-50,1e3*out.q_L,'k-.','LineWidth',lineWidth)
ylabel('Flow rate (L/s)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

leg = legend('p_{p}','p_{L}');
leg.FontSize = fontSize-1;
leg.FontName = 'Cambria Math';

linkaxes([ax1 ax2 ax3],'x')
% xlim([0 600])

sgtitle('Pipeline Boundary Pressure Dynamics',...
    'Interpreter','latex','FontSize',fontSize,'fontname','Times')
% supT = suptitle('WEC-driven Pumping Dynamics');
% supT.FontSize = fontSize+2;
% supT.FontName = 'Times';

end 

%%
% dt = out.par.dt_MOC(1);
% t_vec = 1:floor(par.tend/dt);
% nt = length(t_vec);
% xLP = linspace(0,out.par.L_line,out.par.n_seg(1)+1);
% xHP = linspace(0,out.par.L_line,out.par.n_seg(2)+1);
% 
% f = figure;
% xlabel('Distance from inlet (m)')
% ylabel('Pressure (MPa)')
% 
% for it = 1:nt
%     plot(xLP,1e-6*out.pLP(it,:))
%     ylim([0 2])
%     pause(2*dt)
% end

%%
if 0
xLP = linspace(0,out.par.L_line,out.par.n_seg(1)+1);
xHP = linspace(0,out.par.L_line,out.par.n_seg(2)+1);

figure;
s = surf(xLP,out.t,1e-6*out.pLP);
s.EdgeColor = 'none';
xlabel('Distance from inlet (m)')
ylabel('Time (s)')
zlabel('Pressure (MPa)')


figure;
s = surf(xHP,out.t,1e-6*out.pHP);
s.EdgeColor = 'none';
xlabel('Distance from inlet (m)')
ylabel('Time (s)')
zlabel('Pressure (MPa)')
end
% 
%%
% xLP = linspace(0,out.par.L_line,out.par.n_seg(1)+1);
% xHP = linspace(0,out.par.L_line,out.par.n_seg(2)+1);
% 
% figure;
% s = surf(xLP,out.t,1e3*out.qLP);
% s.EdgeColor = 'none';
% xlabel('Distance from inlet (m)')
% ylabel('Time (s)')
% zlabel('Flow rate (L/s)')
% 
% 
% figure;
% s = surf(xHP,out.t,1e3*out.qHP);
% s.EdgeColor = 'none';
% xlabel('Distance from inlet (m)')
% ylabel('Time (s)')
% zlabel('Flow rate (L/s)')

