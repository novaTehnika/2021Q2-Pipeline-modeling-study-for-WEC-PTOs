function dydt = sys_PTOpipeline(t,y,par,postProcess)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_PTOpipeline.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 3/29/2021
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for a simple wave energy PTO with high 
% and low-pressure pipelines. For use in 2021Q1 pipeline modeling project.
%
% FILE DEPENDENCY: 
% pipeline.m
% nonStateVarsPTOwPL.m
%
% UPDATES:
% 3/29/2021 - created.
% 05/18/2021 - Added unsteady friction model indices to switch
%           statements.
% 05/21/2021 - Revised initalization of dydt for medium line and N Pi-lump 
%           models.
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

%  Purpose: System model; impliments the differential equation for system states
    
    % initalize dy/dt vector array and determine pipeline state indicies
    switch par.pLmodel
        case {1 2 6 7 10 11}
            dydt = zeros(3,1);
        case {3 4 5 8 9}
            dydt = zeros(3 ...
            + ((2 + par.nyf)*par.n_seg(1)-1) ...
            + ((2 + par.nyf)*par.n_seg(2)-1),1);
    end
    
    % extract pressure node states for clarity
    p_hin = y(1);
    p_hout = y(2);
    p_lout = y(3);
    
    % determine non-state, time and state dependent variables
    nonState = nonStateVarsPTOwPL(t,y,par);
    
    % Get states of LP pipeline
    [~, q_lout, dydt_pLineLP,~] = ...
    pipeline(t,y(par.iy_LP),nonState.p_lin,p_lout,par,1,postProcess);  
           
    % Get states of HP pipeline
    [q_hin, q_hout, dydt_pLineHP,~] = ...
    pipeline(t,y(par.iy_HP),p_hin,p_hout,par,2,postProcess); 
    
    % calculate derivitives in pressure states in system nodes
    dp_hin_dt = 1/par.C_hin*(nonState.q_p - q_hin);
    dp_hout_dt = 1/par.C_hout*(q_hout - nonState.q_L);
    dp_lout_dt = 1/par.C_lout*(q_lout - nonState.q_p);
    
    % load dy/dt vector for output
      % system states
    dydt(1) = dp_hin_dt;
    dydt(2) = dp_hout_dt;
    dydt(3) = dp_lout_dt;
    
      % pipeline states (for lumped parameter models)
    switch par.pLmodel
        case {3 4 5 8 9} % long line, lumped parameter model (i.e. n pi-lumps)
            dydt(par.iy_LP) = dydt_pLineLP;
            dydt(par.iy_HP) = dydt_pLineHP;
        otherwise % no pipeline states to load in dy/dt vector
    end
    
end
