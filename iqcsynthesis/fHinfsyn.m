function [K,ga,Xe] = fHinfsyn(Pl,nout,nin,options)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        25-06-2012
% 
% -------------------------------------------------------------------------
%
% Description: This function performs an H-infinity synthesis
%
% Syntax:      [K,gamma]    = fHinfsyn(Pl,nout,nin)
%              [K,gamma]    = fHinfsyn(Pl,nout,nin,options)
%              [K,gamma,Xe] = fHinfsyn(Pl,nout,nin,options)
%
% Usage:       This function solves a particular version of the H-infinity
%              synthesis problem.
%
%              As input:
%
%                # The plant Pl
%
%                              [ A  | Bw1  , Bw2  , Bu  ]
%                       [z1]   [----|-------------------] [w1]
%                  Pl = [z2] = [Cz1 | Dz1w1, Dz1w2, Dz1u] [w2]
%                       [y ]   [Cz2 | Dz2w1, Dz2w2, Dz2u] [u ]
%                              [Cy  | Dyw1 , Dyw2 , Dyu ]
%
%                  with 
%                    - disturbance inputs:  w1, w2
%                    - performance outputs: z1, z2 
%                    - control input:       u
%                    - measurement output:  y
%                    - (A,Bu) stabilizable 
%                    - (A,Cy) detectable
%
%                  Here the first performance channel satisfies:
%
%                  int_0^\infty [*][ I,                  0 ][z1(t)]dt\geq0 
%                               [*][ 0, -options.alpha^2 I ][w1(t)]
%
%                  for all t\geq0 and for for some fixed options.alpha
%                  (default is 1), while the second performance channel
%                  satisfies:
%
%                  int_0^\infty [*][ I,          0 ][z2(t)]dt\geq0
%                               [*][ 0, -gamma^2 I ][w2(t)]
%
%                  for all t\geq0 and where gamma is minimized.
%
%                # The input and output channel data:
%
%                  In the case that nin and nout are defined by
%                    - nin  = nu
%                    - nout = ny
%
%                  it is assumed that nw1 = nz1 = [] and that the remaining
%                  input/output channels are associated with nw2 and nz2
%                  respectively.
%
%                  On the other hand, in the case that nin and nout are
%                  defined by 
%
%                  nin = [nw1;nw2;nu]  and  nout = [nz1;nz2;ny]
%
%                  where:
%
%                    - nw1 denotes the number of disturbance inputs
%                      satisfying performance criteria 1,
%                    - nw2 denotes the number of disturbance inputs
%                      satisfying performance criteria 2,
%                    - nu denotes the number of control inputs.
%
%                  and where
%
%                    - nz1 denotes the number of performance outputs
%                      satisfying performance criteria 1
%                    - nz2 denotes the number of performance outputs
%                      satisfying performance criteria 2
%                    - ny denotes the number of measurements outputs
%
%                # The optional structure "options", where:
%
%                   - options.subopt      If larger than 1, the algorithm
%                                         computes a suboptimal solution
%                                         options.subopt*gamma while
%                                         minimizing the norms on the
%                                         LMI variables (default = 1.01).
%
%                   - options.condnr      If larger than 1, the algorithm
%                                         fixes the suboptimal solution
%                                         options.condnr*gamma while
%                                         minimizing the
%                                         conditioningsnumber of the
%                                         coupling condition [Y,I;I,X] > 0
%                                         (default = 1).
%
%                   - options.uv          Use different options to
%                                         reconstruct the controller. Use
%                                         "options.uv = 1" for U = X and
%                                         V = X^-1-Y and "options.uv = 2"
%                                         for U = Y^-1-X and V = Y
%                                         (default=1). 
%
%                   - options.eliminate   "options.eliminate = 1" computes
%                                         the Hinf controller in one shot
%                                         "options.eliminate = 2" computes
%                                         the Hinf controller by
%                                         eliminating the controller 
%                                         variables first (default = 2).
%
%                   - options.Lyapext     Constructs the extended Lyapunov
%                                         matrix "Xe" in different fashions
%                                         (Default = 5, see fLyapext for
%                                         further details).
%
%                   - options.constants   options.constants = [f1,f2,f3,f4]
%                                         is a verctor whose elements
%                                         perturb the LMIs with LMIi < -fiI
%                                         (Thus fi should be small)
%                                         (default = [0,0,0,0]).
%
%                   - options.init        structure with "options.init.Xe"
%                                         and "options.init.gamma". Specify
%                                         to use the Lyapunov matrix Xe and
%                                         gamma as initial condition in the
%                                         corresponding optimization
%                                         problem. Works only in case
%                                         options.eliminate = 2 and
%                                         options.balreal = 0. Xe should
%                                         have twice the size of A, gamma
%                                         is a scaler (default = []). 
%
%                   - options.Parser      Use LMIlab or Yalmip to solve the
%                                         optimization problem (default =
%                                         'LMIlab') 
%
%                   - options.Solver      Use mincx in case LMIlab is
%                                         considered or another solver
%                                         (e.g. sedumi, sdpt3, mosek, etc.)
%                                         in case Yalmip is used.
%
%                   - options.gmax        Specifies the maximal values of
%                                         gamma Might speed up the algoritm
%                                         (default = 1000)
%
%                   - options.controller  can be used in case to construct
%                                         the controller analitically
%                                         (options.controller = 1, default)
%                                         or by LMI optimization techniques
%                                         (options.controller = 2) 
%
%                   - options.alpha       L2-bound on the first performance
%                                         channel w1->z1 (default = 1).
%
%                   - options.balreal     balance the plant:
%                                          1.) options.balreal = 0: no
%                                              balancing (default)
%                                          2.) options.balreal = 1:
%                                              balreal(sys) (Grammian-based
%                                              IO  balancing)
%                                          3.) options.balreal = 2:
%                                              IQChingbalreal(sys)
%                                              (Hinfinity Grammian-based IO
%                                              balancing. For stable plants
%                                              only.) 
%                                          4.) options.balreal = 3:
%                                              ssbal(sys) (uses "balance"
%                                              to compute a diagonal
%                                              similarity transformation T
%                                              such that [T*A/T , T*B ; C/T
%                                              0] has approximately equal
%                                              row and column norms.
%
%                   - options.FeasbRad    Bound variables (Default = 1e6)
%
%                   - options.Terminate   can be used to change the
%                                         LMIlab-solver details (see help:
%                                         minxc) (default = 20) 
%
%                   - options.RelAcc      can be used to change de
%                                         LMIlab-solver details (see help:
%                                         minxc) (default = 1e-4)
%
%              As Output:
%
%                # The Controller K = ss(Ac,Bc,Cc,Dc) which stabilizes the
%                  closed-loop system lft(Pl,K) while the following
%                  performance criteria are satisfied:
%
%                  int_0^\infty [z1(t)]^T[ I,          0 ][z1(t)]dt\geq0
%                               [w1(t)]  [ 0, -alpha^2 I ][w1(t)]
%
%                  int_0^\infty [z2(t)]^T[ I,       0 ][z2(t)]dt\geq0
%                               [w2(t)]  [ 0, -ga^2 I ][w2(t)]
%
%                  for all t\geq0.
%
%                # The H-infinity norm "ga"
%
%                # The Lyapynov function Xe
%
% -------------------------------------------------------------------------
if nargin == 3
    options.subopt                 = 1.01;
    options.condnr                 = 1;
    options.uv                     = 1;
    options.eliminate              = 2;
    options.Lyapext                = 5;
    options.constants              = [0,0,0,0];
    options.init.Xe                = [];
    options.init.gamma             = [];
    options.Parser                 = 'LMIlab';
    options.Solver                 = 'mincx';
    options.gmax                   = 1000;
    options.FeasbRad               = 1e6;
    options.controller             = 1;
    options.balreal                = 0;
    options.alpha                  = 1;
    options.Terminate              = 20;
    options.RelAcc                 = 1e-4;
elseif nargin == 4
    if isfield(options,'subopt')
        if options.subopt >= 1
            options.subopt         = options.subopt;
        else
            options.subopt         = 1.01;
        end
    else
        options.subopt             = 1.01;
    end
    if isfield(options,'condnr')
        if options.condnr >= 1
            options.condnr         = options.condnr;
        else
            options.condnr         = 1;
        end
    else
        options.condnr             = 1;
    end
    if isfield(options,'alpha')
        if options.alpha>0
            options.alpha          = options.alpha;
        else
            options.alpha          = 1;
        end
    else
        options.alpha              = 1;
    end
    if isfield(options,'uv')
        if options.uv == 1 || options.uv == 2
            options.uv             = options.uv;
        else
            options.uv             = 1;
        end
    else
        options.uv                 = 1;
    end
    if isfield(options,'eliminate')
        if options.eliminate == 1 || options.eliminate == 2
            options.eliminate      = options.eliminate;
        else
            options.eliminate      = 2;
        end
    else
        options.eliminate          = 2;
    end
    if isfield(options,'Lyapext')
        if options.Lyapext == 1 || options.Lyapext == 2 || options.Lyapext == 3 || options.Lyapext == 4 || options.Lyapext == 5
            options.Lyapext        = options.Lyapext;
        else
            options.Lyapext        = 5;
        end
    else
        options.Lyapext            = 5;
    end
    if isfield(options,'constants')
        if min(options.constants) < 0
            options.constants      = [0,0,0,0];
        else
            options.constants      = options.constants;
        end
    else
        options.constants          = [0,0,0,0];
    end
    if isfield(options,'Parser')
        if strcmp(options.Parser,'LMIlab') || strcmp(options.Parser,'Yalmip')
            options.Parser         = options.Parser;
        else
            options.Parser         = 'LMIlab';
        end
    else
        options.Parser             = 'LMIlab';
    end
    if isfield(options,'Solver')
        if strcmp(options.Parser,'LMIlab')
            options.Solver         = 'mincx';
        elseif strcmp(options.Parser,'Yalmip')
            options.Solver         = options.Solver;
        end
    else
        options.Parser             = 'LMIlab';
        options.Solver             = 'mincx';
    end
    if isfield(options,'gmax')
        if options.gmax > 0
            options.gmax           = options.gmax;
        else
            options.gmax           = 1000;
        end
    else
        options.gmax               = 1000;
    end
    if isfield(options,'FeasbRad')
        if options.FeasbRad > 0
            options.FeasbRad       = options.FeasbRad;
        else
            options.FeasbRad       = 1e6;
        end
    else
        options.FeasbRad           = 1e6;
    end
    if isfield(options,'controller')
        if options.controller == 1 || options.controller == 2
            options.controller     = options.controller;
        else
            options.controller     = 1;
        end
    else
        options.controller         = 1;
    end
    if isfield(options,'balreal')
        if options.balreal==0 || options.balreal==1 || options.balreal==2 || options.balreal==3
            options.balreal        = options.balreal;
        else
            options.balreal        = 0;
        end
    else
        options.balreal            = 0;
    end
    if isfield(options,'init')
        if isfield(options.init,'Xe') && isfield(options.init,'gamma') && isfield(options.init,'Tpert')
            if isequal(size(options.init.Xe,1),2*size(Pl.a,1)) && ... 
               isequal(size(options.init.Xe,2),2*size(Pl.a,2)) && ...
               isequal(size(options.init.gamma,1),1) && ...
               isequal(size(options.init.gamma,2),1) && ...
               options.eliminate == 2 && options.solver == 2 && options.balreal == 0
               options.init        = options.init;
            else
                options.init.Xe    = [];
                options.init.Tpert = [];
                options.init.gamma = [];
            end
        elseif isfield(options.init,'Xe') && isfield(options.init,'gamma')
            if isequal(size(options.init.Xe,1),2*size(Pl.a,1)) && ... 
               isequal(size(options.init.Xe,2),2*size(Pl.a,2)) && ...
               isequal(size(options.init.gamma,1),1) && ...
               isequal(size(options.init.gamma,2),1) && ...
               options.eliminate == 2 && options.solver == 2 && options.balreal == 0
               options.init.Tpert          = zeros(size(options.init.Xe,1));
            else
                options.init.Xe    = [];
                options.init.Tpert = [];
                options.init.gamma = [];
            end            
        else
            options.init.Xe        = [];
            options.init.Tpert     = [];
            options.init.gamma     = [];
        end
    else
        options.init.Xe            = [];
        options.init.Tpert         = [];
        options.init.gamma         = [];
    end
    if isfield(options,'Terminate')
        if options.Terminate>0
            options.Terminate      = options.Terminate;
        else
            options.Terminate      = 20;
        end
    else
        options.Terminate          = 20;
    end
    if isfield(options,'RelAcc')
        if options.RelAcc>0
            options.RelAcc         = options.RelAcc;
        else
            options.RelAcc         = 1e-4;
        end
    else
        options.RelAcc             = 1e-4;
    end
end
switch options.balreal
    case 0
        Pl                         = Pl;
    case 1
        Pl                         = balreal(Pl);
    case 2
        Pl                         = fHinfbalreal(Pl);
    case 3
        Pl                         = ssbal(Pl);
end

% Define dimensional info
if size(nin,1) - size(nin,2) == 0
    nw1                            = 0;
    nw2                            = size(Pl.d,2) - nin;
    nu                             = nin;
    nz1                            = 0;
    nz2                            = size(Pl.d,1) - nout;
    ny                             = nout;
else
    if abs(size(Pl.d,2) - sum(nin)) > 0 || abs(size(Pl.d,1) - sum(nout)) > 0
        error('The the vectors "nin" or "nout" are not well-defined');
    end    
    nw1                            = nin(1);
    nw2                            = nin(2);
    nu                             = nin(3);
    nz1                            = nout(1);
    nz2                            = nout(2);
    ny                             = nout(3);
end
nx                                 = size(Pl.a,1);
nw                                 = nw1+nw2;
nz                                 = nz1+nz2;

% Define ss matrices
A                                  = Pl.a;
Bw                                 = Pl.b(:,1:nw);
Bu                                 = Pl.b(:,nw+1:end);
Cz                                 = Pl.c(1:nz,:);
Cy                                 = Pl.c(nz+1:end,:);
Dzw                                = Pl.d(1:nz,1:nw);
Dzu                                = Pl.d(1:nz,nw+1:end);
Dyw                                = Pl.d(nz+1:end,1:nw);
Dyu                                = Pl.d(nz+1:end,nw+1:end);

% Define performance blocks
Pnw1                               = blkdiag(options.alpha*eye(nw1),zeros(nw2));
Pnw2                               = blkdiag(zeros(nw1),eye(nw2));
Pnz1                               = blkdiag(options.alpha*eye(nz1),zeros(nz2));
Pnz2                               = blkdiag(zeros(nz1),eye(nz2));

% Solve the Hinf synthesis problem
switch options.eliminate
    case 1 
        % Case without elimination of the controller variables
        %
        % Define and solve the LMI To'*Tm*To<-he Tr
        %
        %      [0,0,Y,0,0  ,0  ,0,0,0      ,     0 ] [A^T ,0   ,0   ,Cz^T ]
        %      [0,0,0,X,0  ,0  ,0,0,0      ,     0 ] [0   ,I   ,0   ,0    ]
        %      [Y,0,0,0,0  ,0  ,0,0,0      ,     0 ] [I   ,0   ,0   ,0    ]
        %      [0,X,0,0,0  ,0  ,0,0,0      ,     0 ] [0   ,A   ,Bw  ,0    ] 
        %      [0,0,0,0,0  ,0  ,K,L,0      ,     0 ] [0   ,I   ,0   ,0    ]
        % [*]^T[0,0,0,0,0  ,0  ,M,0,0      ,     0 ] [Bu^T,0   ,0   ,Dzu^T] < ...
        %      [0,0,0,0,K^T,M^T,0,0,0      ,     0 ] [I   ,0   ,0   ,0    ]
        %      [0,0,0,0,L^T,0  ,0,0,0      ,     0 ] [0   ,Cy  ,Dyw ,0    ]
        %      [0,0,0,0,0  ,0  ,0,0,-g*Pnw2,     0 ] [0   ,0   ,I   ,0    ]
        %      [0,0,0,0,0  ,0  ,0,0,0      ,-g*Pnz2] [0   ,0   ,0   ,I    ]
        %
        %               [I,0]               [Bu ]               [0,0,0   ,   0]
        %               [0,0][0,A, Bw ,0]   [0  ]               [0,0,0   ,   0]
        %        ...<-he[0,0][0,Cz,Dzw,0]-he[0  ]Dc[0,Cy,Dyw,0]+[0,0,Pnw1,   0]-options.constants(3) I
        %               [0,I]               [Dzu]               [0,0,0   ,Pnw2]
        %
        % together with the coupling constraint
        %
        %   [Y,I] > [options.constants(1) I, 0                     ]
        %   [I,X]   [0                     , options.constants(2) I]

        % Find The direct feedthrough term "Dc":
        if norm(Dzw) == 0
            Dc                     = zeros(nu,ny);
        else
            % Define LMI problem
            probOpt.FeasbRad       = options.FeasbRad;
            probOpt.Parser         = options.Parser;
            probOpt.Solver         = options.Solver;
            probOpt.Terminate      = options.Terminate;
            probOpt.RelAcc         = options.RelAcc;
            probOpt.eps            = 1e-9;
            prob                   = iqcprob(probOpt);
            
            % Define LMI variables
            Dc                     = iqcvar(prob,[nu,ny],'full');
            gamma                  = iqcvar(prob,[1,1],'symmetric');
            
            % Construct Outer factor
            Pnw2var                = blkdiag(zeros(nw1),gamma*eye(nw2));
            Pnz2var                = blkdiag(zeros(nz1),gamma*eye(nz2));
            Tm                     = blkdiag(oblkdiag(Dc),-Pnw2var,-Pnz2var);
            To                     = [blkdiag(Dyw,Dzu');eye(nw+nz)];
            Tc                     = -fOblkdiag(Dzw')+blkdiag(Pnw1,Pnz1);
            
            % Construct LMI
            prob                   = iqclmi(prob,Tm,-1,Tc,To);
            prob                   = iqclmi(prob,gamma,1);
            prob                   = iqcsolve(prob,gamma);
            Dc                     = iqcdec2mat(prob,Dc);
        end
        
        % Define LMI problem
        probOpt.FeasbRad           = options.FeasbRad;
        probOpt.Parser             = options.Parser;
        probOpt.Solver             = options.Solver;
        probOpt.Terminate          = options.Terminate;
        probOpt.RelAcc             = options.RelAcc;
        probOpt.eps                = 1e-9;
        prob                       = iqcprob(probOpt);
        
        % Define LMI variables
        X                          = iqcvar(prob,[nx,nx],'symmetric');
        Y                          = iqcvar(prob,[nx,nx],'symmetric');
        K                          = iqcvar(prob,[nx,nx],'full');
        L                          = iqcvar(prob,[nx,ny],'full');
        M                          = iqcvar(prob,[nu,nx],'full');
        gamma                      = iqcvar(prob,[1,1],'symmetric');
        
        % Construct Outer factor
        To                         = [A',zeros(nx,nx+nw),Cz';fJ(nx,nw+nz,nx)';...
                                     fJ(nx,nw+nz+nx)';zeros(nx,nx),A,Bw,zeros(nx,nz);...
                                     fJ(nx,nw+nz,nx)';Bu',zeros(nu,nx+nw),Dzu';fJ(nx,nw+nz+nx)';...
                                     zeros(ny,nx),Cy,Dyw,zeros(ny,nz);fJt(nw+nz,2*nx)'];
        Tr1                        = -fHe([fJ(nx,nx+nw+nz),fJ(nz,[],2*nx+nw)]*[zeros(nx+nz,nx),[A,Bw;Cz,Dzw],zeros(nx+nz,nz)]);
        Tr2                        = -fHe([Bu;zeros(nx+nw,nu);Dzu]*Dc*[zeros(ny,nx),Cy,Dyw,zeros(ny,nz)]);
        Tr3                        = blkdiag(zeros(2*nx,2*nx),Pnw1,Pnz1);
        Tr                         = Tr1+Tr2+Tr3;
        
        % Construct matrix viaribles
        Tc                         = blkdiag(Y,X);
        Pnw2var                    = blkdiag(zeros(nw1),gamma*eye(nw2));
        Pnz2var                    = blkdiag(zeros(nz1),gamma*eye(nz2));
        Tm                         = blkdiag(oblkdiag(Tc),oblkdiag([K,L;M,zeros(nu,ny)]),-Pnw2var,-Pnz2var);
        
        % Construct LMIs
        prob                       = iqclmi(prob,Tm,-1,Tr,To);
        prob                       = iqclmi(prob,Tc,1,-fOblkdiag(eye(nx)));
        prob                       = iqclmi(prob,gamma,-1,options.gmax);
        
        % solve Hinfsyn problem
        prob                       = iqcsolve(prob,gamma);
        ga                         = prob.gamma;

        if ga == -1
            error('LMIs are not feasible.');
        end

        % Hinf synthesis for sub-optimal gamma
        if options.subopt > 1
            %  Minimize the norm of the variables X,Y,K,L,M subject to the
            %  LMI
            %
            %      [0,0,Y,0,0  ,0  ,0,0] [A^T ,0   ,0   ,Cz^T ]
            %      [0,0,0,X,0  ,0  ,0,0] [0   ,I   ,0   ,0    ]
            %      [Y,0,0,0,0  ,0  ,0,0] [I   ,0   ,0   ,0    ]
            %      [0,X,0,0,0  ,0  ,0,0] [0   ,A   ,Bw  ,0    ] 
            %      [0,0,0,0,0  ,0  ,K,L] [0   ,I   ,0   ,0    ] < ...
            % [*]^T[0,0,0,0,0  ,0  ,M,0] [Bu^T,0   ,0   ,Dzu^T]
            %      [0,0,0,0,K^T,M^T,0,0] [I   ,0   ,0   ,0    ]
            %      [0,0,0,0,L^T,0  ,0,0] [0   ,Cy  ,Dyw ,0    ]
            %
            %                [I,0]                 [Bu ]                 [0,0,0   ,   0]   [0,0,0                  , 0                 ]
            %                [0,0][0,A, Bw ,0]     [0  ]                 [0,0,0   ,   0]   [0,0,0                  , 0                 ]
            %        ...< -he[0,0][0,Cz,Dzw,0] - he[0  ]Dc[0,Cy,Dyw,0] + [0,0,Pnw1,   0] + [0,0,subopt*gampopt*Pnw2, 0                 ] - options.constants(3) I
            %                [0,I]                 [Dzu]                 [0,0,0   ,Pnw2]   [0,0,0                  , subopt*gamopt*Pnz2]
            %
            % and the coupling constraint
            %
            %   [Y,I] > [options.constants(1) I, 0                     ]
            %   [I,X]   [0                     , options.constants(2) I]
        
            % Define LMI problem
            probOpt.FeasbRad       = options.FeasbRad;
            probOpt.Parser         = options.Parser;
            probOpt.Solver         = options.Solver;
            probOpt.Terminate      = options.Terminate;
            probOpt.RelAcc         = options.RelAcc;
            probOpt.eps            = 1e-9;
            prob                   = iqcprob(probOpt);
            
            % Define LMI variables
            X                      = iqcvar(prob,[nx,nx],'symmetric');
            Y                      = iqcvar(prob,[nx,nx],'symmetric');
            K                      = iqcvar(prob,[nx,nx],'full');
            L                      = iqcvar(prob,[nx,ny],'full');
            M                      = iqcvar(prob,[nu,nx],'full');
            gamma                  = iqcvar(prob,[1,1],'symmetric');
            
            % Bound LMIvariables
            prob                   = fBoundvar(prob,X,gamma);
            prob                   = fBoundvar(prob,Y,gamma);
            prob                   = fBoundvar(prob,K,gamma);
            prob                   = fBoundvar(prob,L,gamma);
            prob                   = fBoundvar(prob,M,gamma);
            
            % Construct Outer factors main LMI
            ga                     = options.subopt*ga;
            To                     = [A',zeros(nx,nx+nw),Cz';fJ(nx,nw+nz,nx)';...
                                     fJ(nx,nw+nz+nx)';zeros(nx,nx),A,Bw,zeros(nx,nz);...
                                     fJ(nx,nw+nz,nx)';Bu',zeros(nu,nx+nw),Dzu';fJ(nx,nw+nz+nx)';...
                                     zeros(ny,nx),Cy,Dyw,zeros(ny,nz)];
            Tr1                    = -fHe([fJ(nx,nx+nw+nz),fJ(nz,[],2*nx+nw)]*[zeros(nx+nz,nx),[A,Bw;Cz,Dzw],zeros(nx+nz,nz)]);
            Tr2                    = -fHe([Bu;zeros(nx+nw,nu);Dzu]*Dc*[zeros(ny,nx),Cy,Dyw,zeros(ny,nz)]);
            Tr3                    = blkdiag(zeros(2*nx,2*nx),Pnw1,Pnz1);
            Tr4                    = fJt(nw+nz,2*nx)*ga*blkdiag(Pnw2,Pnz2)*fJt(nw+nz,2*nx)';
            Tr                     = Tr1+Tr2+Tr3+Tr4;
            
            % Construct matrix viaribles
            Tc                     = blkdiag(Y,X);
            Tp                     = oblkdiag([K,L;M,zeros(nu,ny)]);
            Tm                     = blkdiag(oblkdiag(Tc),Tp);
            
            % Define right-hand size LMIs
            Cc                     = -options.constants(3)*eye(size(To,2))+Tr;
            Cm                     = -fOblkdiag((1+options.constants(2))*eye(nx))+blkdiag((options.constants(1))*eye(2*nx));
            
            % Construct LMIs
            prob                   = iqclmi(prob,Tm,-1,Cc,To);
            prob                   = iqclmi(prob,Tc,1,Cm);
           
            % solve Hinfsyn problem
            prob                   = iqcsolve(prob,gamma);
            gambnd                 = prob.gamma;
            
            if gambnd == -1
                 error('Suboptimal (2nd) LMIs are not feasible.');
            end
        end
        if options.condnr > 1
            if options.subopt > 1
                gambnd_subopt      = options.condnr*gambnd;
            else
                nX                 = norm(iqcdec2mat(prob,X));
                nY                 = norm(iqcdec2mat(prob,Y));
                nK                 = norm(iqcdec2mat(prob,K));
                nL                 = norm(iqcdec2mat(prob,L));
                nM                 = norm(iqcdec2mat(prob,M));
                gambnd_subopt      = options.condnr*max([nX,nY,nK,nL,nM]);
                ga                 = options.condnr*ga;
            end
            
            % Define LMI problem
            probOpt.FeasbRad       = options.FeasbRad;
            probOpt.Parser         = options.Parser;
            probOpt.Solver         = options.Solver;
            probOpt.Terminate      = options.Terminate;
            probOpt.RelAcc         = options.RelAcc;
            probOpt.eps            = 1e-9;
            prob                   = iqcprob(probOpt);
            
            % Define LMI variables
            X                      = iqcvar(prob,[nx,nx],'symmetric');
            Y                      = iqcvar(prob,[nx,nx],'symmetric');
            K                      = iqcvar(prob,[nx,nx],'full');
            L                      = iqcvar(prob,[nx,ny],'full');
            M                      = iqcvar(prob,[nu,nx],'full');
            gamma                  = iqcvar(prob,[1,1],'symmetric');
                        
            % Bound LMIvariables
            prob                   = fBoundvar(prob,X,gambnd_subopt);
            prob                   = fBoundvar(prob,Y,gambnd_subopt);
            prob                   = fBoundvar(prob,K,gambnd_subopt);
            prob                   = fBoundvar(prob,L,gambnd_subopt);
            prob                   = fBoundvar(prob,M,gambnd_subopt);
            
            % Construct Outer factor
            To                     = [A',zeros(nx,nx+nw),Cz';fJ(nx,nw+nz,nx)';...
                                     fJ(nx,nw+nz+nx)';zeros(nx,nx),A,Bw,zeros(nx,nz);...
                                     fJ(nx,nw+nz,nx)';Bu',zeros(nu,nx+nw),Dzu';fJ(nx,nw+nz+nx)';...
                                     zeros(ny,nx),Cy,Dyw,zeros(ny,nz)];
            Tr1                    = -fHe([fJ(nx,nx+nw+nz),fJ(nz,[],2*nx+nw)]*[zeros(nx+nz,nx),[A,Bw;Cz,Dzw],zeros(nx+nz,nz)]);
            Tr2                    = -fHe([Bu;zeros(nx+nw,nu);Dzu]*Dc*[zeros(ny,nx),Cy,Dyw,zeros(ny,nz)]);
            Tr3                    = blkdiag(zeros(2*nx,2*nx),Pnw1,Pnz1);
            Tr4                    = fJt(nw+nz,2*nx)*ga*blkdiag(Pnw2,Pnz2)*fJt(nw+nz,2*nx)';
            Tr                     = Tr1+Tr2+Tr3+Tr4;
            
            % Construct matrix viaribles
            Tc                     = blkdiag(Y,X);
            Tp                     = oblkdiag([K,L;M,zeros(nu,ny)]);
            Tm                     = blkdiag(oblkdiag(Tc),Tp);
            Txy                    = [Y,-gamma*eye(nx);-gamma*eye(nx),X];
            
            % Define right-hand size LMIs
            Cc                     = -(options.constants(3))*eye(size(To,2))+Tr;
            Cm                     = -fOblkdiag(eye(nx))+blkdiag(options.constants(1)*eye(nx),options.constants(2)*eye(nx));
            Cxy                    = -fOblkdiag(eye(nx));
            
            % Construct LMIs
            prob                   = iqclmi(prob,Tm,-1,Cc,To);
            prob                   = iqclmi(prob,Tc,1,Cm);
            prob                   = iqclmi(prob,gamma,-1,options.gmax);
            prob                   = iqclmi(prob,Txy,1,Cxy);
            
            % solve Hinfsyn problem
            prob                   = iqcsolve(prob,gamma);
            gam_cond               = prob.gamma;

            if gam_cond == -1
                 error('Suboptimal (3nd) LMIs are not feasible.');
            end
        end

        % call variables
        X                          = iqcdec2mat(prob,X);
        Y                          = iqcdec2mat(prob,Y);
        K                          = iqcdec2mat(prob,K);
        L                          = iqcdec2mat(prob,L);
        M                          = iqcdec2mat(prob,M);
        
        % Construct controller
        switch options.uv
            case 1
                U                  = X;
                V                  = fInv(X)-Y;
            case 2
                U                  = fInv(Y)-X;
                V                  = Y;
        end
        
        % Construct the extended Lyapunov matrix
        Xe                         = fLyapext(X,Y,options.Lyapext);
        Kt                         = [U,X*Bu;zeros(nu,nx),eye(nu)]\([K,L;M,Dc]-blkdiag(X*A*Y,zeros(nu,ny)))/[V',zeros(nx,ny);Cy*Y,eye(ny)];
        Kt(nx+1:end,nx+1:end)      = Dc;
        K                          = (eye(size(Kt,1))+Kt*blkdiag(zeros(nx),Dyu))\Kt;
        K                          = ss(K(1:nx,1:nx),K(1:nx,nx+1:end),K(nx+1:end,1:nx),K(nx+1:end,nx+1:end));
    case 2 
        % Case with elimination of the controller variables
        %
        % First define and solve the LMIs (existence conditions for K)
        % Tap'*Top'*Tmp*Top*Tap < Tap'*Trp*Tap
        % Tad'*Tod'*Tmd*Tod*Tad < Tad'*Trd*Tad
        %
        % which are define as
        %
        %           [0,X, 0     , 0     ][I   ,0    ,0]                [0   ,    0 ,Cz^T ]
        % (x)^T(x)^T[X,0, 0     , 0     ][A   ,Bw   ,0][Psi,0] < -(x)^T[0   , -Pnw1,Dzw^T][Psi,0]-options.constants(3) I
        %           [0,0,-g*Pnw2, 0     ][0   ,I    ,0][  0,I]         [Cz  ,  Dzw ,-Pnz1][  0,I]
        %           [0,0, 0     ,-g*Pnz2][0   ,0    ,I]
        %
        %           [0,Y, 0     , 0     ][-A^T,-Cz^T,0]                [0   ,    0,  Bw ]
        % (x)^T(x)^T[Y,0, 0     , 0     ][I   ,0,    0][Phi,0] >  (x)^T[0   ,-Pnz1,  Dzw][Phi,0]+options.constants(4) I
        %           [0,0, g*Pnz2, 0     ][0   ,I,    0][0  ,I]         [Bw^T,Dzw^T,-Pnw1][  0,I]
        %           [0,0, 0     , g*Pnw2][0   ,0,    I]
        %
        % where Psi=null(Cy,Dyw) and Phi=null(Bu^T,Dzu^T)
        %
        % as well as
        %
        %   [Y,I] > [options.constants(1) I, 0                     ]
        %   [I,X]   [0                     , options.constants(2) I]
        %
        % Then define and solve the LMI (find controller variables)
        %
        % he(diag(cal(X),I,I)*T1)+T4+he(T2*K*T3)<0, where
        %
        %    [A ,0, Bw ,0]     [0, Bu]                      [0,0, 0          , 0          ]
        % T1=[0 ,0, 0  ,0], T2=[I,  0], T3=[0 ,I,  0,0], T4=[0,0, 0          , 0          ]
        %    [0 ,0, 0  ,0]     [0,  0]     [Cy,0,Dyw,0]     [0,0,-Pnw1-g*Pnw2, 0          ]
        %    [Cz,0, Dzw,0]     [0,Dzu]                      [0,0, 0          ,-Pnz1-g*Pnz2]

        % Define LMI problem
        probOpt.FeasbRad           = options.FeasbRad;
        probOpt.Parser             = options.Parser;
        probOpt.Solver             = options.Solver;
        probOpt.Terminate          = options.Terminate;
        probOpt.RelAcc             = options.RelAcc;
        probOpt.eps                = 1e-9;
        prob                       = iqcprob(probOpt);

        % Define LMI variables
        X                          = iqcvar(prob,[nx,nx],'symmetric');
        Y                          = iqcvar(prob,[nx,nx],'symmetric');
        gamma                      = iqcvar(prob,[1,1],'symmetric');
        
        % Construct matrix variables
        Pnw2var                    = blkdiag(zeros(nw1),gamma*eye(nw2));
        Pnz2var                    = blkdiag(zeros(nz1),gamma*eye(nz2));
        Tmp                        = blkdiag(oblkdiag(X),-Pnw2var,-Pnz2var);
        Tmd                        = blkdiag(oblkdiag(Y),Pnz2var,Pnw2var);
        Tc                         = blkdiag(Y,X);
        
        % Construct outer factors
        Top                        = [fJ(nx,nw+nz)';A,Bw,zeros(nx,nz);fJt(nw+nz,nx)'];
        Tod                        = [-A',-Cz',zeros(nx,nw);eye(nx+nw+nz)];
        Trp                        = -fOblkdiag([Cz';Dzw'])+blkdiag(zeros(nx),Pnw1,Pnz1);
        Trd                        = fOblkdiag([Bw;Dzw])-blkdiag(zeros(nx),Pnz1,Pnw1);
        
        % Construct annihilators
        Tap                        = blkdiag(null([Cy,Dyw]),eye(nz));
        Tad                        = blkdiag(null([Bu',Dzu']),eye(nw));
        
        % Perturb LMIs by small delta        
        Ccp                        = 0.5*fHe(Tap'*Trp*Tap);
        Ccd                        = 0.5*fHe(Tad'*Trd*Tad);
        Cm                         = -fOblkdiag(eye(nx));
        Cg                         = options.gmax;
        
        % Construct LMI
        prob                       = iqclmi(prob,Tmp,-1,Ccp,Top*Tap);
        prob                       = iqclmi(prob,Tmd,1,Ccd,Tod*Tad);
        prob                       = iqclmi(prob,Tc,1,Cm);
        prob                       = iqclmi(prob,gamma,-1,Cg);
        
        % Vectorize the initial condition compatible with X and Y
        if isempty(options.init.Xe) && isempty(options.init.gamma)
            prob.Init              = [];
        else
            Xinit                  = options.init.Xe(1:nx,1:nx);
            Tinit                  = chol(fSchur(options.init.Xe-options.init.Tpert,fJ(nx)));
            Tiniti                 = fInv(Tinit);
            Yinit                  = Tiniti*Tiniti';            
            Xvec                   = fVec(Xinit,'sym');
            Yvec                   = fVec(Yinit,'sym');
            prob.Init              = [Xvec;Yvec;options.init.gamma];
        end
        
        % solve Hinfsyn problem
        prob                       = iqcsolve(prob,gamma);
        ga                         = prob.gamma;
        
        if ga == -1
            error('LMIs are not feasible.');
        end

        if options.subopt > 1
            %                                                                                                 [0   ,    0 ,Cz^T ]
            % (x)^T(x)^T[0,X][I   ,0    ,0][Psi,0] < -(x)^T(x)^T[-g*Pnw2, 0     ][0   ,I    ,0][Psi,0] - (x)^T[0   , -Pnw1,Dzw^T][Psi,0]-options.constants(3) I                
            %           [X,0][A   ,Bw   ,0][  0,I]              [ 0     ,-g*Pnz2][0   ,0    ,I][  0,I]        [Cz  ,  Dzw ,-Pnz1][  0,I]
            %                                              
            %                                                                                              [0   ,    0,  Bw ]
            % (x)^T(x)^T[0,Y][-A^T,-Cz^T,0][Phi,0] > -(x)^T(x)^T[g*Pnz2, 0     ][0   ,I,    0][Phi,0] (x)^T[0   ,-Pnz1,  Dzw][Phi,0]+options.constants(4) I
            %           [Y,0][I   ,0,    0][0  ,I]              [0     , g*Pnw2][0   ,0,    I][0  ,I]      [Bw^T,Dzw^T,-Pnw1][  0,I]
            %            
            % where Psi=null(Cy,Dyw) and Phi=null(Bu^T,Dzu^T)
            
            % Define LMI problem
            probOpt.FeasbRad       = options.FeasbRad;
            probOpt.Parser         = options.Parser;
            probOpt.Solver         = options.Solver;
            probOpt.Terminate      = options.Terminate;
            probOpt.RelAcc         = options.RelAcc;
            probOpt.eps            = 1e-9;
            prob                   = iqcprob(probOpt);
            
            % Define LMI variables
            X                      = iqcvar(prob,[nx,nx],'symmetric');
            Y                      = iqcvar(prob,[nx,nx],'symmetric');
            gamma                  = iqcvar(prob,[1,1],'symmetric');
                        
            % Bound LMIvariables
            prob                   = fBoundvar(prob,X,gamma);
            prob                   = fBoundvar(prob,Y,gamma);
            
            % Construct matrix variables
            Tmp                    = oblkdiag(X);
            Tmd                    = oblkdiag(Y);
            Tc                     = blkdiag(Y,X);
            
            % Construct outer factors
            ga                     = options.subopt*ga;
            Top                    = [fJ(nx,nw+nz)';A,Bw,zeros(nx,nz)];
            Tod                    = [-A',-Cz',zeros(nx,nw);fJ(nx,nw+nz)'];
           
            Trp                    = ga*fJt(nw+nz,nx)*blkdiag(Pnw2,Pnz2)*fJt(nw+nz,nx)'-fOblkdiag([Cz';Dzw'])+blkdiag(zeros(nx),Pnw1,Pnz1);
            Trd                    = -ga*[zeros(nw+nz,nx),eye(nw+nz)]'*blkdiag(Pnz2,Pnw2)*[zeros(nw+nz,nx),eye(nw+nz)]+fOblkdiag([Bw;Dzw])-blkdiag(zeros(nx),Pnz1,Pnw1);
            
            % Construct annihilators
            Tap                    = blkdiag(null([Cy,Dyw]),eye(nz));
            Tad                    = blkdiag(null([Bu',Dzu']),eye(nw));
            
            % Perturb LMIs by small delta
            Ccp                    = -options.constants(3)*eye(size(Tap,2))+0.5*fHe(Tap'*Trp*Tap);
            Ccd                    = options.constants(4)*eye(size(Tad,2))+0.5*fHe(Tad'*Trd*Tad);
            Cm                     = -fOblkdiag(eye(nx))+blkdiag((options.constants(1))*eye(nx),(options.constants(2))*eye(nx));
            
            % Construct LMI
            prob                   = iqclmi(prob,Tmp,-1,Ccp,Top*Tap);
            prob                   = iqclmi(prob,Tmd,1,Ccd,Tod*Tad);
            prob                   = iqclmi(prob,Tc,1,Cm);
            
            % solve Hinfsyn problem
            prob                   = iqcsolve(prob,gamma);
            gambnd                 = prob.gamma;
            
            if gambnd == -1
                 error('Suboptimal (2nd) LMIs are not feasible!');
            end
        end
        if options.condnr > 1
            if options.subopt > 1
                gambnd_subopt      = options.condnr*gambnd;
            else
                nX                 = norm(iqcdec2mat(prob,X));
                nY                 = norm(iqcdec2mat(prob,Y));
                gambnd_subopt      = options.condnr*max([nX,nY])^2;
                ga                 = options.condnr*ga;
            end
            
            % Define LMI problem
            probOpt.FeasbRad       = options.FeasbRad;
            probOpt.Parser         = options.Parser;
            probOpt.Solver         = options.Solver;
            probOpt.Terminate      = options.Terminate;
            probOpt.RelAcc         = options.RelAcc;
            probOpt.eps            = 1e-9;
            prob                   = iqcprob(probOpt);
                       
            % Define LMI variables
            X                      = iqcvar(prob,[nx,nx],'symmetric');
            Y                      = iqcvar(prob,[nx,nx],'symmetric');
            gamma                  = iqcvar(prob,[1,1],'symmetric');
            
            % Bound LMIvariables
            prob                   = fBoundvar(prob,X,gambnd_subopt);
            prob                   = fBoundvar(prob,Y,gambnd_subopt);
            
            % Construct matrix variables
            Tmp                    = oblkdiag(X);
            Tmd                    = oblkdiag(Y);
            Tc                     = blkdiag(Y,X);
            Txy                    = [Y,-gamma*eye(nx);-gamma*eye(nx),X];
            
            % Construct outer factors
            Top                    = [fJ(nx,nw+nz)';A,Bw,zeros(nx,nz)];
            Tod                    = [-A',-Cz',zeros(nx,nw);fJ(nx,nw+nz)'];
            Trp                    = -fOblkdiag([Cz';Dzw'])+blkdiag(zeros(nx),Pnw1+ga*Pnw2,Pnz1+ga*Pnz2);
            Trd                    = fOblkdiag([Bw;Dzw])-blkdiag(zeros(nx),Pnz1+ga*Pnz2,Pnw1+ga*Pnw2);
            
            % Construct annihilators
            Tap                    = blkdiag(null([Cy,Dyw]),eye(nz));
            Tad                    = blkdiag(null([Bu',Dzu']),eye(nw));
            
            % Perturb LMIs by small delta
            Ccp                    = -options.constants(3)*eye(size(Tap,2))+0.5*fHe(Tap'*Trp*Tap);
            Ccd                    = options.constants(4)*eye(size(Tad,2))+0.5*fHe(Tad'*Trd*Tad);
            Cm                     = -fOblkdiag(eye(nx))+blkdiag(options.constants(1)*eye(nx),options.constants(2)*eye(nx));
            Cxy                    = -fOblkdiag(eye(nx));
            
            % Construct LMI
            prob                   = iqclmi(prob,Tmp,-1,Ccp,Top*Tap);
            prob                   = iqclmi(prob,Tmd,1,Ccd,Tod*Tad);
            prob                   = iqclmi(prob,Txy,1,Cm);
            prob                   = iqclmi(prob,Txy,1,Cxy);
            
            % solve Hinfsyn problem
            prob                   = iqcsolve(prob,gamma);
            gam_cond               = prob.gamma;

            if gam_cond == -1
                 error('Suboptimal (3nd) LMIs are not feasible.');
            end
        end
        
        % call variables
        X                          = iqcdec2mat(prob,X);
        Y                          = iqcdec2mat(prob,Y);
        
        % Construct the extended Lyapunov matrix
        Xe                         = fLyapext(X,Y,options.Lyapext);
        % disp(['||Xe||: ', num2str(norm(Xe))]);
                
        % Construct the synthesis LMIs and substitute Xe and gamma
        Xee                        = blkdiag(Xe,eye(nw+nz));
        T1                         = Xee*[A,zeros(nx),Bw,zeros(nx,nz);zeros(nx+nw,2*nx+nw+nz);Cz,zeros(nz,nx),Dzw,zeros(nz)];
        T2                         = Xee*[fJ(nx,nw+nz,nx),[Bu;zeros(nx+nw,nu);Dzu]];
        T3                         = [fJ(nx,nw+nz,nx)';Cy,zeros(ny,nx),Dyw,zeros(ny,nz)];
        T4                         = blkdiag(zeros(2*nx),-Pnw1-ga*Pnw2,-Pnz1-ga*Pnz2);
        switch options.controller
            case 1 % Analytically construct the controller
                Kt                 = fInvproj(T2',T3,fHe(T1)+T4);
            case 2 % Find the controller by optimization
                if norm(Dzw) == 0
                    Dc             = zeros(nu,ny);
                else
                    % Define LMI problem
                    probOpt.FeasbRad  = options.FeasbRad;
                    probOpt.Parser    = options.Parser;
                    probOpt.Solver    = options.Solver;
                    probOpt.Terminate = options.Terminate;
                    probOpt.RelAcc    = options.RelAcc;
                    probOpt.eps       = 1e-9;
                    prob              = iqcprob(probOpt);
                    
                    % Define LMI variable
                    Dc             = iqcvar(prob,[nu,ny],'full');
                    gamma          = iqcvar(prob,[1,1],'symmetric');

                    % Construct Outer factor
                    Pnw2var        = blkdiag(zeros(nw1),gamma*eye(nw2));
                    Pnz2var        = blkdiag(zeros(nz1),gamma*eye(nz2));
                    Tm             = blkdiag(oblkdiag(Dc),-Pnw2var,-Pnz2var);
                    To             = [blkdiag(Dyw,Dzu');eye(nw2+nz2)];
                    Cc             = -fOblkdiag(Dzw')+blkdiag(Pnw1,Pnz1);
                    
                    % Construct LMI
                    prob           = iqclmi(prob,Tm,-1,Cc,To);
                    prob           = iqclmi(prob,gamma,1);
                    
                    % solve optimization problem
                    prob           = iqcsolve(prob,gamma);
                    Dc             = iqcdec2mat(prob,Dc);
                end
                NT2                = null(T2');
                NT3                = null(T3);
                
                % Check feasibility of the first LMIs
                Ep                 = eig(NT2'*(fHe(T1)+T4)*NT2);
                Ed                 = eig(NT3'*(fHe(T1)+T4)*NT3);
                if max(Ep) > 0 || max(Ed) > 0
                    disp(Ep);
                    disp(Ed);
                    error('The main analysis LMI is not satisfied')
                end
                
                % Create LMI problem
                probOpt.FeasbRad   = options.FeasbRad;
                probOpt.Parser     = options.Parser;
                probOpt.Solver     = options.Solver;
                probOpt.Terminate  = options.Terminate;
                probOpt.RelAcc     = options.RelAcc;
                probOpt.eps        = 1e-9;
                prob               = iqcprob(probOpt);
                
                % Define LMI variables
                K                  = iqcvar(prob,[nx,nx],'full');
                L                  = iqcvar(prob,[nx,ny],'full');
                M                  = iqcvar(prob,[nu,nx],'full');
                gamma              = iqcvar(prob,[1,1],'symmetric');

                % Bound LMIvariables
                prob               = fBoundvar(prob,K,gamma);
                prob               = fBoundvar(prob,L,gamma);
                prob               = fBoundvar(prob,M,gamma);
                
                % Create LMI
                P                  = oblkdiag([K,L;M,zeros(nu,ny)]);
                Co                 = -fHe(T1)-fHe(T2*[zeros(nx,nx+ny);zeros(nu,nx),Dc]*T3)-T4;
                prob               = iqclmi(prob,P,-1,Co,[T2';T3]);
                
                % Solve LMI problem
                prob               = iqcsolve(prob,gamma);
                K                  = iqcdec2mat(prob,K);
                L                  = iqcdec2mat(prob,L);
                M                  = iqcdec2mat(prob,M);
                
                % Construct the Controller
                Kt                 = [K,L;M,Dc];
        end
K                                  = (eye(size(Kt,1))+Kt*blkdiag(zeros(nx),Dyu))\Kt;
K                                  = ss(K(1:nx,1:nx),K(1:nx,nx+1:end),K(nx+1:end,1:nx),K(nx+1:end,nx+1:end));
D                                  = fCleanupsys(ss(K.d),1e-14);
K.d                                = D.d;
end
end