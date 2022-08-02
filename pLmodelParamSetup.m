function par = pLmodelParamSetup(par,n_seg)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pLmodelParamSetup.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 3/25/2021
%
% PURPOSE/DESCRIPTION:
% Setup variables used in and in dealing with pipeline model in the
% required format.
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 3/25/2021 - created.
% 4/16/2021 - Changed determination of number of segments and MOC 
% time step to be based on the number of segments (n_seg) instead of time 
% step (dt_moc).
% 05/18/2021 - Added unsteady friction model indices to switch
%           statements.
% 05/21/2021 - Updated model number assignments. Revised state indices 
%           (par.iy_LP and par.iy_HP) for medium line and N Pi-lump models.
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
        case {1 2} % short line
            par.n_seg(1) = 1;
            par.n_seg(2) = 1;
            par.iy_LP = [];
            par.iy_HP = [];

        case 3 % medium line spoofing short line
            par.n_seg(1) = 1;
            par.n_seg(2) = 1;
            par.iy_LP = 3 + (1:2*par.n_seg(1)-1);
            par.iy_HP = 3 + (2*par.n_seg(1):...
                        (2*par.n_seg(1)-1) + (2*par.n_seg(2)-1));

            par.I(1) = par.rho*(0.1)/par.A_line;
            par.I(2) = par.rho*(0.1)/par.A_line;
            
        case {4 8} % medium line
            par.n_seg(1) = 1;
            par.n_seg(2) = 1;
            par.iy_LP = 3 + (1:(2 + par.nyf)-1);
            par.iy_HP = 3 + ((2 + par.nyf):((2 + par.nyf)-1) + ((2 + par.nyf)-1));

            par.I(1) = par.rho*(par.L_line)/par.A_line;
            par.I(2) = par.rho*(par.L_line)/par.A_line;
            
        case {5 9} % N Pi-lump
            par.n_seg(1) = n_seg;
            par.n_seg(2) = n_seg; 
            par.iy_LP = 3 + (1:(2 + par.nyf)*par.n_seg(1)-1);
            par.iy_HP = 3 + ((2 + par.nyf)*par.n_seg(1): ...
                        ((2 + par.nyf)*par.n_seg(1)-1) + ...
                        ((2 + par.nyf)*par.n_seg(2)-1));


            par.I(1) = par.rho*(par.L_line/par.n_seg(1))/par.A_line;
            par.I(2) = par.rho*(par.L_line/par.n_seg(2))/par.A_line;
            
        case {6 7 10 11}  % MOC

            % number of segments for each line
%             par.n_seg(1) = floor(par.L_line/(par.dt_MOC(1)*par.a(1))); % LP line
%             par.n_seg(2) = floor(par.L_line/(par.dt_MOC(2)*par.a(2))); % HP line

            par.n_seg(1) = n_seg; % LP line
            par.n_seg(2) = n_seg; % HP line
            
             % MOC solver time step (total timestep)
             par.dt_MOC(1) = par.L_line/par.n_seg(1)/par.a(1); % [s] LP line
             par.dt_MOC(2) = par.L_line/par.n_seg(2)/par.a(2); % [s] HP line


            % ideal gas law constant
            for ID = 1:2
                par.C1(ID) = (par.p_o-par.p_vap)*par.R*...
                              par.A_line*par.L_line/(par.n_seg(ID)-1);
            end

            par.iy_LP = [];
            par.iy_HP = [];

    end

end