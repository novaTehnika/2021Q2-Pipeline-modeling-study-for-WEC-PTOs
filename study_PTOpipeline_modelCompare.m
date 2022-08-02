% study_PTOpipeline_modelCompare.m script m-file
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

rng(2); % set the seed for the random number generator
    
%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% specify simualtion parameters
designCase = 4;
designSubCase = 3;

pLmodel = [3 4 5 6 7];
nVar = length(pLmodel);

for i = 1:nVar

    disp([num2str(i),' of ',num2str(nVar)]) % output to console 
                                            % for user monitoring
                                            
    par = parameters_PTOpipeline(par,designCase,designSubCase,pLmodel(i));

    % run simulation
    tic
    out(i) = sim_PTOpipeline(par,1);
    toc

end
%% %%%%%%%%%%%%   POST PROCESS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract data
for i = 1:nVar
    
    PPlossLP(i) = 1e-3*out(i).PPmean_LP_fric;
    PPlossHP(i) = 1e-3*out(i).PPmean_HP_fric;
    
    p_loutMu(i) = 1e-6*out(i).p_loutDist.mu;
    p_hinMu(i) = 1e-6*out(i).p_hinDist.mu;
    p_houtinMu(i) = 1e-6*out(i).p_houtDist.mu;
    
    p_loutSigma(i) = 1e-6*out(i).p_loutDist.sigma;
    p_hinSigma(i) = 1e-6*out(i).p_hinDist.sigma;
    p_houtinSigma(i) = 1e-6*out(i).p_houtDist.sigma;
    
    j = find(out(i).dp_loutdtDist.f > .997,1);
    dp_loutdt99x7(i) = 1e-6*out(i).dp_loutdtDist.xi(j);
    j = find(out(i).dp_hindtDist.f > .997,1);
    dp_hindt99x7(i) = 1e-6*out(i).dp_hindtDist.xi(j);
    j = find(out(i).dp_houtdtDist.f > .997,1);
    dp_houtdt99x7(i) = 1e-6*out(i).dp_houtdtDist.xi(j);
    
    deltap_wecMu(i) = 1e-6*out(i).deltap_wecDist.mu;
    deltap_wecSigma(i) = 1e-6*out(i).deltap_wecDist.sigma;
end

% calcualte error wrt to the DGCM
error = @(x,i) 100*abs(x(i) - x(end))/x(end);

for i = 1:nVar
    
    e_PPlossLP(i) = error(PPlossLP,i);
    e_PPlossHP(i) = error(PPlossHP,i);
    
    e_p_loutMu(i) = error(p_loutMu,i);
    e_p_hinMu(i) = error(p_hinMu,i);
    e_p_houtinMu(i) = error(p_houtinMu,i);
    
    e_p_loutSigma(i) = error(p_loutSigma,i);
    e_p_hinSigma(i) = error(p_hinSigma,i);
    e_p_houtinSigma(i) = error(p_houtinSigma,i);
    
    e_deltap_wecMu(i) = error(deltap_wecMu,i);
    e_deltap_wecSigma(i) = error(deltap_wecSigma,i);
    

    e_dp_loutdt99x7(i) = error(dp_loutdt99x7,i);
    e_dp_hindt99x7(i) = error(dp_hindt99x7,i);
    e_dp_houtdt99x7(i) = error(dp_houtdt99x7,i);
    
end

precVal = 3;
precError = 2;
for i = 1:nVar
    
    A(i,1) = PPlossLP(i);
    A(i,2) = PPlossHP(i);
    A(i,3) = p_loutMu(i);
    A(i,4) = p_hinMu(i);
    A(i,5) = p_houtinMu(i);
    A(i,6) = p_loutSigma(i);
    A(i,7) = p_hinSigma(i);
    A(i,8) = p_houtinSigma(i);
    A(i,9) = dp_loutdt99x7(i);
    A(i,10) = dp_hindt99x7(i);
    A(i,11) = dp_houtdt99x7(i);
    A(i,12) = deltap_wecMu(i);
    A(i,13) = deltap_wecSigma(i);
    
end
A = A';

for i = 1:nVar
    
    B(i,1) = e_PPlossLP(i);
    B(i,2) = e_PPlossHP(i);
    B(i,3) = e_p_loutMu(i);
    B(i,4) = e_p_hinMu(i);
    B(i,5) = e_p_houtinMu(i);
    B(i,6) = e_p_loutSigma(i);
    B(i,7) = e_p_hinSigma(i);
    B(i,8) = e_p_houtinSigma(i);
    B(i,9) = e_dp_loutdt99x7(i);
    B(i,10) = e_dp_hindt99x7(i);
    B(i,11) = e_dp_houtdt99x7(i);
    B(i,12) = e_deltap_wecMu(i);
    B(i,13) = e_deltap_wecSigma(i);
    
end
B = B';

for i = 1:nVar
    
    C(i,1) = {[num2str(PPlossLP(i),precVal),' (',num2str(e_PPlossLP(i),precError),')']};
    C(i,2) = {[num2str(PPlossHP(i),precVal),' (',num2str(e_PPlossHP(i),precError),')']};
    C(i,3) = {[num2str(p_loutMu(i),precVal),' (',num2str(e_p_loutMu(i),precError),')']};
    C(i,4) = {[num2str(p_hinMu(i),precVal),' (',num2str(e_p_hinMu(i),precError),')']};
    C(i,5) = {[num2str(p_houtinMu(i),precVal),' (',num2str(e_p_houtinMu(i),precError),')']};
    C(i,6) = {[num2str(p_loutSigma(i),precVal),' (',num2str(e_p_loutSigma(i),precError),')']};
    C(i,7) = {[num2str(p_hinSigma(i),precVal),' (',num2str(e_p_hinSigma(i),precError),')']};
    C(i,8) = {[num2str(p_houtinSigma(i),precVal),' (',num2str(e_p_houtinSigma(i),precError),')']};
    C(i,9) = {[num2str(dp_loutdt99x7(i),precVal),' (',num2str(e_dp_loutdt99x7(i),precError),')']};
    C(i,10) = {[num2str(dp_hindt99x7(i),precVal),' (',num2str(e_dp_hindt99x7(i),precError),')']};
    C(i,11) = {[num2str(dp_houtdt99x7(i),precVal),' (',num2str(e_dp_houtdt99x7(i),precError),')']};
    C(i,12) = {[num2str(deltap_wecMu(i),precVal),' (',num2str(e_deltap_wecMu(i),precError),')']};
    C(i,13) = {[num2str(deltap_wecSigma(i),precVal),' (',num2str(e_deltap_wecSigma(i),precError),')']};
    
end
C = C';
%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PPlossLPmean(i) = out(i).PPmean_LP_fric;
% PPlossLPmean(i) = out(i).PPmean_HP_fric;
% 
% out.p_loutDist = out.p_loutDist
% out.p_hinDist
% out.p_houtDist
% 
% out.dp_loutdtDist
% out.dp_hindtDist
% out.dp_houtdtDist
% 
% out.deltap_wecDist
% 
% out.pLP = pLP;
% out.pHP = pHP;
