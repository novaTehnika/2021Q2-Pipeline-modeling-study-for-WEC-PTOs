% study_PTOpipeline_modelCompare_uFric.m script m-file
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
% 4/17/2021 - Added post-processing for design metrics.
% 05/18/2021 - Updated model versions implimenting unsteady 
%           friction in clear cmd (2nd executable line).
% 10/18/2021 - Moved rng(*) function to the inside of the simulation loop
% so that the random number generator seed is the same for every simulation
% in the study.
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
par.tend = 1200; %[s] time span of simulation

par.odeSolverRelTol = 1e-6; % Rel. error tolerance parameter for ODE solver
par.odeSolverAbsTol = 1e-6; % Abs. error tolerance parameter for ODE solver


    
%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% specify simualtion parameters
designCase = 6;
designSubCase = 1;

pLmodel = [5 9 6 10];
nVar1 = length(pLmodel);
nVar2 = 29;
switch designCase
    case 5 % switching at the pump
        Tsw = logspace(log10(0.05),log10(10),nVar2);
        L_line = 1000;
        fl = 1466/(2*L_line);
        fp = 2./Tsw;
        fp/fl
    case 6 % switching at the load
        Tsw = logspace(log10(0.05),log10(10),nVar2);
        L_line = 1000;
        fl = 1466/(2*L_line);
        fp = 2./Tsw;
        fp/fl
end
% nVar2 = length(Tsw);

for i = 1:nVar1
    parfor j = 1:nVar2
        disp(['Run (',num2str(i),',',num2str(j),') of ('...
               ,num2str(nVar1),',',num2str(nVar2),')']) % output to console 
                                                      % for user monitoring
                                                      
        rng(2); % set the seed for the random number generator
        parSim = parameters_PTOpipelineStudy(par,designCase,designSubCase,pLmodel(i));
        
        switch designCase
            case 5 % switching at the pump
                parSim.Tsw_pump = Tsw(j);    
            case 6 % switching at the load
                parSim.Tsw_load = Tsw(j);    
        end

        % run simulation
        tic
        out(i,j) = sim_PTOpipeline(parSim,1);
        toc
        
    end
end
%% %%%%%%%%%%%%   POST PROCESS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract data
for i = 1:nVar1
    for j = 1:nVar2
        
    PPlossLP(i,j) = 1e-3*out(i,j).PPmean_LP_fric;
    PPlossLP_sfric(i,j) = 1e-3*out(i,j).PPmean_LP_sfric;
    PPlossLP_ufric(i,j) = 1e-3*out(i,j).PPmean_LP_ufric;
    
    PPlossHP(i,j) = 1e-3*out(i,j).PPmean_HP_fric;
    PPlossHP_sfric(i,j) = 1e-3*out(i,j).PPmean_HP_sfric;
    PPlossHP_ufric(i,j) = 1e-3*out(i,j).PPmean_HP_ufric;
    
    p_loutMu(i,j) = 1e-6*out(i,j).p_loutDist.mu;
    p_hinMu(i,j) = 1e-6*out(i,j).p_hinDist.mu;
    p_houtinMu(i,j) = 1e-6*out(i,j).p_houtDist.mu;
    
    p_loutSigma(i,j) = 1e-6*out(i,j).p_loutDist.sigma;
    p_hinSigma(i,j) = 1e-6*out(i,j).p_hinDist.sigma;
    p_houtinSigma(i,j) = 1e-6*out(i,j).p_houtDist.sigma;
    
    k = find(out(i,j).dp_loutdtDist.f > .997,1);
    dp_loutdt99x7(i,j) = 1e-6*out(i,j).dp_loutdtDist.xi(k);
    k = find(out(i,j).dp_hindtDist.f > .997,1);
    dp_hindt99x7(i,j) = 1e-6*out(i,j).dp_hindtDist.xi(k);
    k = find(out(i,j).dp_houtdtDist.f > .997,1);
    dp_houtdt99x7(i,j) = 1e-6*out(i,j).dp_houtdtDist.xi(k);
    
    deltap_wecMu(i,j) = 1e-6*out(i,j).deltap_wecDist.mu;
    deltap_wecSigma(i,j) = 1e-6*out(i,j).deltap_wecDist.sigma;
    
    end
end

% calcualte error wrt to the DGCM
error = @(x,i) 100*abs(x(i,j) - x(end))/x(end);

for i = 1:nVar1
    for j = 1:nVar2
    
    e_PPlossLP(i,j) = error(PPlossLP,i);
    e_PPlossHP(i,j) = error(PPlossHP,i);
    
    e_p_loutMu(i,j) = error(p_loutMu,i);
    e_p_hinMu(i,j) = error(p_hinMu,i);
    e_p_houtinMu(i,j) = error(p_houtinMu,i);
    
    e_p_loutSigma(i,j) = error(p_loutSigma,i);
    e_p_hinSigma(i,j) = error(p_hinSigma,i);
    e_p_houtinSigma(i,j) = error(p_houtinSigma,i);
    
    e_deltap_wecMu(i,j) = error(deltap_wecMu,i);
    e_deltap_wecSigma(i,j) = error(deltap_wecSigma,i);
    

    e_dp_loutdt99x7(i,j) = error(dp_loutdt99x7,i);
    e_dp_hindt99x7(i,j) = error(dp_hindt99x7,i);
    e_dp_houtdt99x7(i,j) = error(dp_houtdt99x7,i);
    
    end
end


%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% total power loss in high-pressure pipeline
figure

semilogx(Tsw,PPlossHP(1,:),'kx','MarkerSize',10,'LineWidth',2)
hold on
semilogx(Tsw,PPlossHP(2,:),'rx','MarkerSize',10,'LineWidth',2)

semilogx(Tsw,PPlossHP(3,:),'ko','MarkerSize',10,'LineWidth',2)
semilogx(Tsw,PPlossHP(4,:),'ro','MarkerSize',10,'LineWidth',2)

grid on 

yLim = ylim;
yLim(1) = 0;
ylim(yLim)

legend('N Pi-lump','N Pi-lump with unsteady fric.','fMOC','fMOC with unsteady fric.')
ylabel('Power loss (kW)')
xlabel('Switching period (s)')
title('Total Flow Loss: High-Pressure Pipeline')

%% unsteady friction power loss in high-pressure pipeline
figure

semilogx(Tsw,PPlossHP_ufric(2,:),'o','MarkerSize',10,'LineWidth',2)
hold on
semilogx(Tsw,PPlossHP_ufric(4,:),'o','MarkerSize',10,'LineWidth',2)

grid on 

yLim = ylim;
yLim(1) = 0;
ylim(yLim)

legend('N Pi-lump with unsteady fric.','fMOC with unsteady fric.')
ylabel('Power loss (kW)')
xlabel('Switching period (s)')
title('Unsteady Flow Loss: High-Pressure Pipeline')

%% unsteady friction power loss in high-pressure pipeline as a function 
% of the ratio between swithcing and line freq.
L_line = 1000;
fl = 1466/(2*L_line);
fp = 2./Tsw;
fp/fl

figure

semilogx(fp/fl,PPlossHP_ufric(2,:),'o','MarkerSize',10,'LineWidth',2)
hold on
semilogx(fp/fl,PPlossHP_ufric(4,:),'o','MarkerSize',10,'LineWidth',2)

grid on 

yLim = ylim;
yLim(1) = 0;
ylim(yLim)

ylabel('Power loss (kW)')
xlabel('Swithcing freq./Line freq. (Hz/Hz)')
title('Unsteady Flow Loss: High-Pressure Pipeline')

%% unsteady friction power loss in high-pressure pipeline
figure


semilogx(Tsw,PPlossHP_sfric(2,:),'kx','MarkerSize',10,'LineWidth',2)
hold on
semilogx(Tsw,PPlossHP_sfric(4,:),'ko','MarkerSize',10,'LineWidth',2)

semilogx(Tsw,PPlossHP_ufric(2,:),'rx','MarkerSize',10,'LineWidth',2)
hold on
semilogx(Tsw,PPlossHP_ufric(4,:),'ro','MarkerSize',10,'LineWidth',2)

grid on 

yLim = ylim;
yLim(1) = 0;
ylim(yLim)


legend('Steady, N Pi-lump with unsteady fric.','Steady, fMOC with unsteady fric.',...
    'Unsteady, N Pi-lump with unsteady fric.','Unsteady, fMOC with unsteady fric.')
ylabel('Power loss (kW)')
xlabel('Switching period (s)')
title('Flow Loss: High-Pressure Pipeline')


%% Womersley number

% Womersley_wave = par.d_line/2*(2*par.Tp.^-1*par.rho/par.mu).^0.5
% Womersley_sw = par.d_line/2*(Tsw.^-1*par.rho/par.mu).^0.5

figure
semilogx(Tsw,Womersley_sw)
grid on
ylabel('Womersley number')
xlabel('Switching period (s)')
title('Womersley Number Versus Switching Freqency')

