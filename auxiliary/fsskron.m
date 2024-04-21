function F = fsskron(G,nr)
% -------------------------------------------------------------------------
%
% IQClab:      Version 3.4.0
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available for non-commercial usage under a
%              Creative Commons (Attribution-NoDerivatives 4.0
%              International (CC BY-ND 4.0)) license:  
%              https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
% Author:      J.Veenman
% Date:        09-12-2019
% 
% -------------------------------------------------------------------------
%
% Description: This function computes the Kronecker product
%
%                            [kron(G.a,Inr) | kron(G.b,Inr)]
%              kron(G,Inr) = [--------------|--------------]
%                            [kron(G.c,Inr) | kron(G.d,Inr)]
% 
%              with
%              1.) G is state space object
%              2.) nr is the number of repetitions of G
%
% Syntax:      F = fsskron(G,nr)
%
% Usage:       F = fsskron(G,nr) creates a new state space object F given
%              the input G and the number of repetitions nr.
%
% -------------------------------------------------------------------------

Inr = eye(nr);
A   = kron(G.a,Inr);
B   = kron(G.b,Inr);
C   = kron(G.c,Inr);
D   = kron(G.d,Inr);
F   = ss(A,B,C,D,G.Ts);
end