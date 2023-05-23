function [q_in, q_out, dydt_pLine, pLsolnOut] = ...
    	 pipelineMOC(t,y,p_in,p_out,par,lineID,postProcess)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pipelineMOC.m function m-file
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
% FILE DEPENDENCY:
% mocSimple.m
% mocDGCM.m
% mocSimple_uFric.m
% mocDGCM_uFric.m
%
% UPDATES:
% 3/25/2021 - created.
% 5/17/2021 - added logic and functions for MOC models with
%           unsteady friction. Note that the unsteady model requires 2dt
%           worth of information, longer arrays are set up for these
%           models.
% 05/21/2021 - Updated model number assignments.
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

    if t == 0 && ~ postProcess; clear pL pLsolnLocal; end
    
    persistent pL pLsolnLocal %#ok<*TLEV>
    
    if t == 0 && ~ postProcess % initialize pL and pLsolnLocal structures
        % for all pipelines
        for id = 1:par.n_lines
            pL(id).t_MOC = 0;   % time for MOC solutions
            pL(id).nt = 1;      % number of timesteps set up
                                % (length of arrays)
            pL(id).itLast_MOC = pL(id).nt; % memory of last MOC
                                                    % iteration solved
            pLsolnLocal(id, 1).q = []; % stored MOC solutions
            pLsolnLocal(id, 1).p = [];
            pLsolnLocal(id, 1).PPfric = 0;
            if par.pLmodel == 7 || par.pLmodel == 11
                pLsolnLocal(id, 1).qd = [];
                pLsolnLocal(id, 1).qu = [];
                pLsolnLocal(id, 1).V = [];
            end
            if par.pLmodel == 10 || par.pLmodel == 11
                if par.pLmodel == 10
                    pLsolnLocal(id, 1).y = [];
                elseif par.pLmodel == 11
                    pLsolnLocal(id, 1).yu = [];
                    pLsolnLocal(id, 1).yd = [];
                end
                pLsolnLocal(id, 1).PPs = 0;
                pLsolnLocal(id, 1).PPu = 0;
            end
        end
    end
    
    if pL(lineID).nt == 1 % set initial conditions
        it = 1;
        switch par.pLmodel
            case {6 7}
                nX = par.n_seg(lineID) + 1; % number of grid points saved 
                                            % outside of the MOC function
            case {10 11}
                nX = 2*par.n_seg(lineID) + 1;
        end
        pLsolnLocal(lineID, it).q = zeros(nX, 1);
        p_ave = (p_in + p_out) / 2;
        pLsolnLocal(lineID, it).p = ...
          p_ave * ones(nX, 1);
      
        if par.pLmodel == 7  || par.pLmodel == 11
            V = par.C1(lineID) / (p_ave - par.p_vap);
            pLsolnLocal(lineID, it).qd = ...
              zeros(2 * par.n_seg(lineID), 1);
            pLsolnLocal(lineID, it).qu = ...
              zeros(2 * par.n_seg(lineID), 1);
            pLsolnLocal(lineID, it).V = ...
              V * ones(2 * par.n_seg(lineID) - 1, 1);
        end
        if par.pLmodel == 10 || par.pLmodel == 11
            if par.pLmodel == 10
                pLsolnLocal(lineID, 1).y = ...
                    zeros(2*par.n_seg(lineID) + 1, length(par.n));
            elseif par.pLmodel == 11
                pLsolnLocal(lineID, 1).yu = ...
                    zeros(2*par.n_seg(lineID), length(par.n));
                pLsolnLocal(lineID, 1).yd = ...
                    zeros(2*par.n_seg(lineID), length(par.n));
            end
        end
        
    end

	switch postProcess % determine function call mode
                       % (0 - solution mode, 1 - post-process)

        case 0 % solution mode: determine if update is need
               % and what states to output

            % if: ODE solver has passed threshold for MOC, MOC needs to 
            %     update for new MOC time step
            if t >= pL(lineID).t_MOC(pL(lineID).itLast_MOC) ...
                 + par.dt_MOC(lineID)

                % update MOC time
                pL(lineID).t_MOC(pL(lineID).nt + 1) ...
                 = pL(lineID).t_MOC(end) + par.dt_MOC(lineID);
                it = pL(lineID).nt + 1;
                pL(lineID).nt = it;

                % set flag to solve MOC
                solveFlag = 1;

                % elseif: the ODE solver has gone back in time.
                % Output older MOC solution.
            elseif t < pL(lineID).t_MOC(pL(lineID).itLast_MOC)

                % determine what time step for the MOC solution
                % is appropriate
                it = find((pL(lineID).t_MOC(:) < t), 1, 'last');

                % set flag to solve MOC
                solveFlag = 0;

                % else: ODE solver is between MOC solutions.
                % No need to solve. Output results of last solve.
            else
                it = pL(lineID).itLast_MOC;
                solveFlag = 0;

            end % end if (time step checking)

            % if: solver flag is raised, solve if needed. Solve and update 
            %     memory.
            if solveFlag

                %  update last call memory
                pL(lineID).itLast_MOC = it;

                % Call MOC function for state update
                switch par.pLmodel
                    case 6
                        pLsolnLocal(lineID, it) = ...
                            mocSimple(t, ...
                            pLsolnLocal(lineID, it - 1), ...
                            p_in, p_out, lineID, par);
                    case 7
                        pLsolnLocal(lineID, it) = ...
                            mocDGCM(t, ...
                            pLsolnLocal(lineID, it - 1), ...
                            p_in, p_out, lineID, par);
                    case 10
                        pLsolnLocal(lineID, it) = ...
                            mocSimple_uFric(t, ...
                            pLsolnLocal(lineID, it - 1), ...
                            p_in, p_out, lineID, par);
                    case 11
                        pLsolnLocal(lineID, it) = ...
                            mocDGCM_uFric(t, ...
                            pLsolnLocal(lineID, it - 1), ...
                            p_in, p_out, lineID, par);
                end

            end % end if (solution routine)

        case 1 % post-process mode: determine MOC time index for stored
           %                    results

            % determine index for the stored MOC solution
            it = find(pL(lineID).t_MOC(:) <= t, 1, 'last');

	end % end switch (solution versus post-process modes)

    % output selected solution results
    q_in = pLsolnLocal(lineID, it).q(1);
    q_out = pLsolnLocal(lineID, it).q(end);
    dydt_pLine = [];
    switch par.pLmodel
        case {6 7}
            pLsolnOut.p = pLsolnLocal(lineID, it).p(:);
            pLsolnOut.q = pLsolnLocal(lineID, it).q(:);
        case {10 11}
            pLsolnOut.p = pLsolnLocal(lineID, it).p(1:2:end);
            pLsolnOut.q = pLsolnLocal(lineID, it).q(1:2:end);
            pLsolnOut.PPs = pLsolnLocal(lineID, it).PPs;
            pLsolnOut.PPu = pLsolnLocal(lineID, it).PPu;
    end
    pLsolnOut.PPfric = pLsolnLocal(lineID, it).PPfric;
    
end
   