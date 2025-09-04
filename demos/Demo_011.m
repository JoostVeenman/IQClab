% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        17-01-2020
% 
% -------------------------------------------------------------------------
% Demo_011:    This demo demonstrates how to:
%
%                1.) Create an LMI problem
%                2.) Implement LMI constraints
%                3.) Solve an LMI problem
%                4.) Perform a solution check
% -------------------------------------------------------------------------
clc;close all;clear ga

% Specify inputs
clc
disp('------------------------------------------------------------------');
disp('Demo_011: Implementation of LMIs:');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo demonstrates how to:');
disp(' ');
disp('  1.) Create an LMI problem');
disp('  2.) Implement LMI constraints');
disp('  3.) Solve an LMI problem');
disp('  4.) Perform a solution check');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('Size of the problem - Please specify the dimension of the matrix T');
disp(' ');
n = input('Choose an integer n [nr. LMI variables is approx. 0.5*n^2]: ');
disp(' ');
disp('------------------------------------------------------------------');

% Define random LTI system with n states
G       = rss(n,1,2);
[no,ni] = size(G);

% Define IQC-LMI problem
prob  = iqcprob; %('Parser','Yalmip','Solver','sdpt3');

% Define LMI variables
X     = iqcvar(prob,[n,n],'symmetric');
gamma = iqcvar(prob,[1,1],'full');

% Define 1st LMI
P     = blkdiag(oblkdiag(X),-gamma*eye(no+ni));
A     = [eye(n),zeros(n,no+ni);G.a,G.b,zeros(n,no);zeros(ni+no,n),eye(no+ni)];
B     = fHe([zeros(n+ni,n+no+ni);-[G.c,G.d,zeros(no)]]);

prob  = iqclmi(prob,P,-1,B,A);
prob  = iqclmi(prob,X,1);

% Solve IQC-LMI problem
prob  = iqcsolve(prob,gamma);

% Check solution
P     = iqcdec2mat(prob,P);
X     = iqcdec2mat(prob,X);
gamma = iqcdec2mat(prob,gamma);

disp('Check if the main LMI is negative definite:');
disp(eig(A'*P*A-B));
disp('Check if X is positive definite');
disp(eig(X));
disp('The computed bound on the H-infinity norm is:');
disp(gamma);



