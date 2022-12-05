% Author: Lucas de Sousa Rosa
% Date: 05/12/2022
%
% ------------------------------------------------------------------------
% 
% Attempt to find a basic feasible solution for the linear program of the form
% 
%                           minimize c'x 
%                         subject to Ax = b
%                                    x >= 0
%
% solving the auxiliary problem. The auxiliary problem is of the form
%
%                           minimize 0'x + 1'y
%                         subject to Ax = b 
%                                    x >= 0
%                                    y >= 0
%
% ------------------------------------------------------------------------
%
%   Input parameters:
%
% A - (m,n) constraint matrix 
% b - (m,1) right-hand side vector
% c - (n,1) cost vector
% m - number of constraints
% n - number of variables
%
%   Output parameters:
%
% ind - status parameter
%   ind =  0 indicates a basic feasible solution was found
%   ind =  1 the problem is infeasible
%  
% x - (m,1) basic feasible solution found
% bind - (m,1) indices of basic variables
% nbind - (n-m,1) indices of nonbasic variables
%

function [ind,x,bind,nbind] = simplex_init(A,b,c,m,n)
    
endfunction