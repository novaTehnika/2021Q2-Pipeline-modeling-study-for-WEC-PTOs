function [q_in, q_out, dydt_pLine, pLsoln] = ...
         pipelineRlam(t,y,p_in,p_out,par,lineID,postProcess)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pipelineRlam.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 3/25/2021
%
% PURPOSE/DESCRIPTION:
% 
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 3/25/2021 - created.
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

    % resistance assuming laminar flow
    R = 8*pi*par.L_line/par.A_line^2; 
    
    % solve for flow
    q = (p_in - p_out)/R; 

    % model outputs
    q_in = q;
    q_out = q;
    dydt_pLine = [];

    % Calculate friction loss
    if postProcess
        pLsoln.PPfric = R*q^2; 
    else
        pLsoln = [];
    end   
    
end