function iqcprob = iqcsolve(iqcprob,minvar)

% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        26-11-2019
% 
% -------------------------------------------------------------------------
%
% Description: solve lmi problem
%
% Syntax:      iqcprob = iqcsolve(iqcprob)
%              iqcprob = iqcsolve(iqcprob,minvar)
%
% Usage:       iqcprob = iqcsolve(iqcprob,minvar) minimizes the LMI
%              variable minvar subject to the LMI problem iqcprob.
%
%              Similarly, iqcprob = iqcsolve(iqcprob) check the feasibility
%              of the LMI problem iqcprob.
%
% -------------------------------------------------------------------------
if iqcprob.BoundVars < 2
    if strcmp(iqcprob.Parser,'Yalmip')     % solve IQC analysis problem with Yalmip
        if strcmp(iqcprob.Display,'off') || strcmp(iqcprob.Display,'offe')
            Displ = 0;
        elseif strcmp(iqcprob.Display,'on')
            Displ = 1;
        end
        Options = sdpsettings('solver',iqcprob.Solver,'verbose',Displ,'radius',iqcprob.FeasbRad);
        if nargin == 1
            if ~strcmp(iqcprob.Display,'offe')
                disp('-------------------------------------------------------------------');
                text1 = ['LMI feasibility problem with ' num2str(yalmip('nvars')) ' descision variables.'];
                disp(text1);
                disp(' - Parser: Yalmip');
                text2 = [' - Solver: ' iqcprob.Solver];
                disp(text2);
                disp('-------------------------------------------------------------------');
            end
            diagnostics = optimize(iqcprob.lmi,[],Options);
            if strcmp(iqcprob.Display,'on')
                disp('-------------------------------------------------------------------');
            end
            if diagnostics.problem == 0
                iqcprob.gamma = 'feasible';
            else
                iqcprob.gamma = 'infeasible';
            end
        elseif nargin == 2
            if ~strcmp(iqcprob.Display,'offe')
                disp('-------------------------------------------------------------------');
                text1 = ['LMI optimization problem with ' num2str(yalmip('nvars')) ' descision variables.'];
                disp(text1);
                disp(' - Parser: Yalmip');
                text2 = [' - Solver: ' iqcprob.Solver];
                disp(text2);
                disp('-------------------------------------------------------------------');
            end
            optimize(iqcprob.lmi,minvar.var,Options);
            if strcmp(iqcprob.Display,'on')
                disp('-------------------------------------------------------------------');
            end
            if double(minvar.var) == 0
                iqcprob.gamma = -1;
            else
                iqcprob.gamma = double(minvar.var);
            end
        end
    elseif strcmp(iqcprob.Parser,'LMIlab') % solve IQC analysis problem with LMIlab
        if strcmp(iqcprob.Display,'off') || strcmp(iqcprob.Display,'offe')
            Displ = 1;
        elseif strcmp(iqcprob.Display,'on')
            Displ = 0;
        end
        Options = [iqcprob.RelAcc,iqcprob.MaxNumIter,iqcprob.FeasbRad,iqcprob.Terminate,Displ];
        if nargin == 1
            iqcprob.lmi = getlmis;
            if ~strcmp(iqcprob.Display,'offe')
                disp('-------------------------------------------------------------------');
                text = ['LMI feasibility problem with ' num2str(decnbr(iqcprob.lmi)) ' descision variables.'];
                disp(text);
                disp(' - Parser: LMIlab');
                disp(' - Solver: feasp');
                disp('-------------------------------------------------------------------');
            end
            [tmin,xsol] = feasp(iqcprob.lmi,Options);
            if strcmp(iqcprob.Display,'on')
                disp('-------------------------------------------------------------------');
            end
            if tmin > -1e-4
                iqcprob.gamma = 'infeasible';
                iqcprob.xsol = [];
            else
                iqcprob.gamma = 'feasible';
                iqcprob.xsol = xsol;
            end
        elseif nargin == 2
            iqcprob.lmi = getlmis;
            c = zeros(decnbr(iqcprob.lmi),1);
            if ~strcmp(iqcprob.Display,'offe')
                disp('-------------------------------------------------------------------');
                text = ['LMI optimization problem with ' num2str(decnbr(iqcprob.lmi)) ' descision variables.'];
                disp(text);
                disp(' - Parser: LMIlab');
                disp(' - Solver: mincx');
                disp('-------------------------------------------------------------------');
            end
            c(minvar.nvar) = 1;
            if isempty(iqcprob.Init)
                [gamma,xsol] = mincx(iqcprob.lmi,c,Options,[],0);
            else
                if ~strcmp(iqcprob.Display,'offe')
                    disp('Solve LMI problem with initial condition');
                end
                [gamma,xsol] = mincx(iqcprob.lmisys,c,iqcprob.options,iqcprob.Init);
            end
            if strcmp(iqcprob.Display,'on')
                disp('-------------------------------------------------------------------');
            end
            if isempty(gamma) || isempty(xsol)
                iqcprob.gamma = -1;
                iqcprob.xsol = [];
            else
                iqcprob.gamma = gamma;
                iqcprob.xsol = xsol;
            end            
        end
    end    
end
end
