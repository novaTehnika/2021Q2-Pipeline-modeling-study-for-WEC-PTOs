function nonState = nonStateVarsPTOwPL(t,y,par)
% nonStateVarsPTOwPL.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 05/03/2021
%
% PURPOSE/DESCRIPTION:
% Calculate non-state variable based on time and system state
%
% FILE DEPENDENCY:
% pumpFlow.m
%
% UPDATES:
% 05/03/2021 - Created.
% 06/01/2021 - Added swithching logic to pump and load flow 
%            rates.
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

    % extract pressure node states for clarity
    % p_hin = y(1); % [Pa] pressure at outlet of HP pipeline
    p_hout = y(2); % [Pa] pressure at inlet of HP pipeline
    % p_lout = y(3); % [Pa] pressure at outlet of LP pipeline

    % pressure at the inlet to the LP pipeline
    nonState.p_lin = par.p_tank; 

    % time-dependent flow rate from driven by the pump
    nonState.q_p = pumpFlow(t,par) ...
        *(mod(t,par.Tsw_pump)/par.Tsw_pump < 1 - par.duty_pump); 

    % flow rate through the system load
    nonState.q_L = p_hout/par.R_load ...
        *(mod(t,par.Tsw_load)/par.Tsw_load < 1 - par.duty_load); 

end
