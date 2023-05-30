function [q_in, q_out, dydt_pLine, pLsoln] = ...
         pipeline(t,y,p_in,p_out,par,lineID,postProcess)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pipeline.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 03/29/2021
%
% PURPOSE/DESCRIPTION:
% A switch for pipeline models based on the parameter par.PLmodel.
%
% FILE DEPENDENCY:
% pipelineRlam.m
% pipelineRturb.m
% pipelinePi.m
% pipelineNPi.m
% pipelinePi_uFric.m
% pipelineNPi_uFric.m
% pipelineMOC.m
%     mocSimple.m
%     mocDGCM.m
%     mocSimple_uFric.m
%     mocDGCM_uFric.m
%
% UPDATES:
% 03/29/2021 - created.
% 05/21/2021 - Updated model number assignments to include 
%            models with unsteady friction. Added functions for lumped 
%            parameter models with unsteady friction.
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

    switch par.pLmodel
        case 1 % short line with laminar friction
            [q_in, q_out, dydt_pLine, pLsoln] = ...
            pipelineRlam(t,y,p_in,p_out,par,lineID,postProcess);

        case 2 % short line with laminar and turbulent friction
            [q_in, q_out, dydt_pLine, pLsoln] = ...
            pipelineRturb(t,y,p_in,p_out,par,lineID,postProcess);

        case {3 4} % medium line model
            [q_in, q_out, dydt_pLine, pLsoln] = ...
            pipelinePi(t,y,p_in,p_out,par,lineID,postProcess);

        case 5 % long line model (n pi-lumps)
            [q_in, q_out, dydt_pLine, pLsoln] = ...
            pipelineNPi(t,y,p_in,p_out,par,lineID,postProcess);
        
        case 8 % medium line model with unsteady friction
            [q_in, q_out, dydt_pLine, pLsoln] = ...
            pipelinePi_uFric(t,y,p_in,p_out,par,lineID,postProcess);
        
        case 9 % long line model (n pi-lumps) with unsteady
            [q_in, q_out, dydt_pLine, pLsoln] = ...
            pipelineNPi_uFric(t,y,p_in,p_out,par,lineID,postProcess);

        case {6 7 10 11} % MOC models with and without unsteady friction
            [q_in, q_out, dydt_pLine, pLsoln] = ...
            pipelineMOC(t,y,p_in,p_out,par,lineID,postProcess);
        
    end

end
