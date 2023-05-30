% study_PTOpipeline_MOCDGCMConvergence.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 03/29/2021
%
% PURPOSE/DESCRIPTION:
% This script executes simulations for a specified design case (a 
% selection of system parameters specified in parameters_PTOpipelinStudy.m)  
% using the fixed wavelength, constant parameter Method of Characteristics
% pipeline model augmented with the descrete gas cavity model. The
% grid spacing is varied for the purpose of testing convergence.
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
designCase = 2
designSubCase = 4
pLmodel = 6
par = parameters_PTOpipelineStudy(par,designCase,designSubCase,pLmodel);


% study parameters
% dt_moc = (5e0)*0.5.^(0:1:8);
% nVar = length(dt_moc);
n_seg = (2)*2.^(0:1:7);
nVar = length(n_seg);

% initialize convergence metrics
PPbalLP = zeros(nVar,1);
PPbalHP = zeros(nVar,1);
PPmean_LP_fric = zeros(nVar,1);
PPmean_HP_fric = zeros(nVar,1);
VbalLP = zeros(nVar,1);
VbalHP = zeros(nVar,1);
dp_houtdt_0x997pTile = zeros(nVar,1);
deltap_wecMean = zeros(nVar,1);
deltap_wecStd = zeros(nVar,1);


for i = 1:nVar
    
    disp([num2str(i),' of ',num2str(nVar)]) % output to console 
                                            % for user monitoring
                                            
    par = pLmodelParamSetup(par,n_seg(i));

   tic
    out(i) = sim_PTOpipeline(par,0);
   toc

end

%%
for i = 1:nVar
    
    PPbalLP(i) = out(i).PPbalLP;
    PPbalHP(i) =  out(i).PPbalHP;
    PPmean_LP_fric(i) = out(i).PPmean_LP_fric;
    PPmean_HP_fric(i) = out(i).PPmean_HP_fric;

    VbalLP(i) = out(i).VbalLP;
    VbalHP(i) = out(i).VbalHP;

    iPtile = find(out(i).dp_houtdtDist.f < 0.997,1,'last');
    dp_houtdt_0x997pTile(i) = out(i).dp_houtdtDist.xi(iPtile);
    deltap_wecMean(i) = out(i).deltap_wecDist.mu;
    deltap_wecStd(i) = out(i).deltap_wecDist.sigma;

end

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%
% Power and Power Balance Error
%%%%%%
figure
yyaxis left
loglog(n_seg,1e-3*PPmean_LP_fric)
hold on
loglog(n_seg,1e-3*abs(PPbalLP-PPmean_LP_fric))
ylabel('Power (kW)')
yyaxis right
loglog(n_seg,100*abs(PPbalLP./PPmean_LP_fric))
ylabel('Error x100%')
legend('Model','Boundary')
title('Power Loss, LP')
xlabel('Number of pipeline segments')


figure
yyaxis left
loglog(n_seg,1e-3*PPmean_HP_fric)
hold on
loglog(n_seg,1e-3*abs(PPbalHP-PPmean_HP_fric))
ylabel('Power (kW)')
yyaxis right
loglog(n_seg,100*abs(PPbalLP./PPmean_HP_fric))
ylabel('Error x100%')
legend('Model','Boundary')
title('Power Loss, HP')
xlabel('Number of pipeline segments')

%%%%%%
% Volume Balance
%%%%%%
figure
semilogx(n_seg,VbalLP)
hold on
semilogx(n_seg,VbalHP)
title('Volume Balance')
xlabel('Number of pipeline segments')
ylabel('Volume (m^3)')


%%%%%%%
% dp_houtdt_0x997pTile
%%%%%%
figure

semilogx(n_seg,1e-6*dp_houtdt_0x997pTile)
title('Rate of Change in Pressure of p_{h,out}')
xlabel('Number of pipeline segments')
ylabel('Rate of Change in Pressure (MPa/s)')


%%%%%%
% deltap_wecMean
%%%%%%
figure
semilogx(n_seg,1e-6*deltap_wecMean)
title('Mean Pressure Differentional Across WEC-Driven Pump')
xlabel('Number of pipeline segments')
ylabel('Pressure (MPa)')



% title(['Design Case = ',num2str(designCase),'.',num2str(designSubCase),', PL Model = ',num2str(pLmodel)])
