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
  
  prinf("Fase 1\n")
  
  % Guarantees viability
  for i = 1:m
    if b(i) < 0
      A(i,:) = -A(i,:);
      b(i) = -b(i);
    endif
  endfor
  
  % Introduce auxiliary variables
  D = [A eye(m)];
  c_aux = [zeros(n,1); ones(m,1)];
  bind = n+1:n+m;
  nbind = 1:n;
  x = zeros(n,1);
  x(bind) = b;
  
  ind = 2; % ensure we execute at least one step
  
  % Simplex phase two iterations
  while ind == 2
    [ind,x,bind,nbind,Binv] = simplex_step(D,b,c_aux,m,n+m,x,bind,nbind,Binv)
  endwhile
  
  % Check for positive optimal cost in the auxiliary problem
  if c_aux(bind)' * x(bind) > 0:
    ind = 1;
    return
  else
    art = intersect(bind,n+1:n+m)
    while ~isempty(art) % checks for artifical variables in the basis
      l = art(1);
      us = Binv * A;
      if´
    endwhile
    ind == 0;
    nbind = setdiff(1:n,bind);
    x(nbind) = zeros(n-m,1);
    return
  endif
  
endfunction