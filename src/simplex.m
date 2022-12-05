% Author: Lucas de Sousa Rosa
% Date: 05/12/2022
%
% ------------------------------------------------------------------------
% 
% Two phase revised simplex method for solving linear programming problems of
% the form
% 
%                           minimize c'x 
%                         subject to Ax = b
%                                    x >= 0
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
%   ind = -1 the problem is unbounded
%   ind =  0 the problem has optimal solution
%   ind =  1 the problem is infeasible
%  
% x - (n,1) vector specifying the solution
%   x is the last basic feasible solution found (when feasible)
%
% d - (m,1) vector specifying the direction
%   d is the direction on which the objective function goes to -infinity (when
%   unbounded)
%
function [ind, x, d] = simplex(A,b,c,m,n) 
  
endfunction