% Author: Lucas de Sousa Rosa
% Date: 05/12/2022
%
% ------------------------------------------------------------------------
% 
% Performs a single step of the simplex method for the linear program of the 
% form
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
% x - (m,1) current basic feasible solution
% bind - (m,1) indices of current basic variables
% nbind - (n-m,1) indices of current nonbasic variables
% Binv - (m,m) inverse matrix of the current basis B
%
%   Output parameters:
%
% ind - status parameter
%   ind = -1 the problem is unbounded
%   ind =  0 the problem has optimal solution
%   ind =  2 simplex method step completed
%  
% x - (m,1) basic feasible solution after one step
% bind - (m,1) indices of basic variables after one step
% nbind - (n-m,1) indices of nonbasic variables after one step
% Binv - (m,m) inverse matrix of the basis B after one step
%
function [ind,x,bind,nbind,Binv] = simplex_step(A,b,c,m,n,x,bind,nbind,Binv)
  printf('Variáveis básicas:\n');
  for i = 1:m
    printf('x%d = %f\n',bind(i),x(bind(i)));
  endfor
  printf('\n');
  
  printf('Custo atual da função objetivo: %f\n\n',c'*x);
  
  % Compute vector of simplex multipliers (z = p')
  z = c(bind)'*Binv;
  
  % Compute reduced costs for non-basic variables (basic variables have reduced cost 0)
  rc = c(nbind) - (z*A(:,nbind))';
  
  printf('Custos reduzidos (não-básicas):\n');
  for i = 1:length(nbind)
    printf('rc(x%d) = %f\n',nbind(i),rc(i));
  endfor
  printf('\n');
  
  % Check optimality
  if isempty(rc(rc < 0))
      ind = 0;
      return
  endif
  
  % Select entering variable using Bland's rule
  for i = 1:length(nbind)
    if rc(i) < 0
      j = nbind(i);
      break
    endif
  endfor
    
  printf('Entra na base: x%d\n\n',j);
  
  % Compute direction vector
  u = Binv*A(:,j);
  
  % Check unboundedness
  if all(u <= 0)
      ind = -1;
      % Calculate unbounded direction
      d = zeros(n,1);
      d(j) = 1;
      d(bind) = -u;
      x = d;
      return
  endif
  
  printf('Vetor de direção (índices básicos):\n');
  for i = 1:m
    printf('d%d = %f\n',bind(i),-u(i));
  endfor
  printf('\n');
  
  % Compute step length (select leaving variable using Bland's rule)
  l = 0;
  theta_star = inf;
  for i = 1:m
    if u(i) > 0
      theta = x(bind(i))/u(i);
      if theta < theta_star || (abs(theta - theta_star) < eps && bind(i) < bind(l))
        theta_star = theta;
        l = i;
      endif
    endif
  endfor
  
  printf('Theta*: %f\n\n', theta_star);
  printf('Sai da base: x%d\n\n',bind(l));
  
  % Update basic solution
  x(j) = theta_star;
  x(bind) = x(bind) - theta_star*u;

  % Update basic and non-basic indexes
  bind(l) = j;
  nbind = setdiff(1:n,bind);
  
  % Update inverse of the basis matrix using elementary row operations
  Binv_u = [Binv u];

  % Elementary row operations
  Q = eye(m);
  Q(:, l) = -u ./ u(l);
  Q(l, l) = 1 / u(l);
  Binv_u = Q * Binv_u;
  
  % Update inverse of the basis matrix
  Binv = Binv_u(:, 1:m);
  
  % simplex method step completed
  ind = 2;
  return
endfunction

%!test
%! A = [1 2 2 1 0 0; 2 1 2 0 1 0; 2 2 1 0 0 1];
%! b = [20; 20; 20];
%! c = [-10; -12; -12; 0; 0; 0];
%! m = 3;
%! n = 6;
%! x = [0; 0; 0; 20; 20; 20];
%! bind = [4; 5; 6];
%! nbind = setdiff(1:n,bind);
%! Binv = inv([A(:,bind)]);
%! [ind,x,bind,nbind,Binv] = simplex_step(A,b,c,m,n,x,bind,nbind,Binv);
%! assert(ind, 2);
%! assert(x, [10; 0; 0; 10; 0; 0]);
%! assert(bind, [4; 1; 6]);
%! assert(nbind, setdiff(1:n,[4; 1; 6]))
%! assert(Binv, [1.0 -0.5 0; 0 0.5 0; 0 -1.0 1])