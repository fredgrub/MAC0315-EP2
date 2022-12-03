function [ind v] = simplex(A,b,c,m,n,x,bind,Binv)
% Revised Simplex method (phase two) for solving the linear programming problem
% min c'x subject to Ax = b, x >= 0
%
% Input: A, b, c, m, n, x, bind, Binv
% A = m x n constraint matrix
% b = m x 1 right-hand side vector
% c = n x 1 cost vector
% m = number of constraints
% n = number of variables
% x = n x 1 initial basic feasible solution
% bind = m x 1 vector of indices of basic variables
% Binv = m x m inverse of the basis matrix
%
% Output: ind, v
% ind = 0 if the solution is optimal
% ind = -1 if the solution is unbounded
% v = optimal solution if ind = 0
% v = unbounded direction if ind = -1

printf('Método Simplex Revisado (fase 2)\n')
iter = 0;
rule = "min";

while 1
  printf('*** Iteração #%d ***\n\n',iter);
  iter = iter + 1;

  % Display the indices of the basic variables and the respective values of the basic variables
  printf('Variáveis básicas:\n');
  for i = 1:m
    printf('x%d = %f\n',bind(i),x(bind(i)));
  end
  printf('\n');

  % Display current objective function value
  printf('Custo atual da função objetivo: %f\n\n',c'*x);

  % Compute vector of simplex multipliers (z = p')
  z = c(bind)'*Binv;

  % Compute non-basic index set
  nonbind = setdiff(1:n,bind);

  % Compute reduced costs for non-basic variables (basic variables have reduced cost 0)
  rc = c(nonbind) - (z*A(:,nonbind))';
  
  % Display reduced costs and indices of non-basic variables
  printf('Custos reduzidos (não-básica):\n');
  for i = 1:length(nonbind)
    printf('rc(x%d) = %f\n',nonbind(i),rc(i));
  end
  printf('\n');

  % Check optimality
  if isempty(rc(rc < 0))
      ind = 0;
      v = x;
      printf('Solução ótima encontrada. Valor da função objetivo: %f\n',c'*x);
      for i = 1:n
        printf('x%d = %f\n',i,x(i));
      end
      return
  end

  % Select entering variable (using Bland's rule or minimum reduced cost rule)
  if strcmp(rule, "bland") == 1
    for i = 1:length(nonbind)
      if rc(i) < 0
        j = nonbind(i);
        break
      end
    end
  elseif strcmp(rule, "min") == 1
    j = nonbind(min(find(rc == min(rc))));
  endif
 
  % Display entering variable
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
      v = d;
      printf('Solução Ilimitada.\n');
      printf('Direção ilimitada:\n');
      for i = 1:n
        printf('d%d = %f\n',i,d(i));
      end
      return
  end

  % Display the indices of the basic variables and the respective values of the direction vector
  printf('Vetor de direção (apenas índices básicos):\n');
  for i = 1:m
    printf('d%d = %f\n',bind(i),-u(i));
  end
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
      end
    end
  end

  % Display step length
  printf('Theta*: %f\n\n', theta_star);

  % Display leaving variable
  printf('Sai da base: x%d\n\n',bind(l));

  % Update basic solution
  x(j) = theta_star;
  x(bind) = x(bind) - theta_star*u;

  % Update basis
  bind(l) = j;

  % Update inverse of the basis matrix using elementary row operations
  Binv_u = [Binv u];

  % Elementary row operations
  Q = eye(m);
  Q(:, l) = -u ./ u(l);
  Q(l, l) = 1 / u(l);
  Binv_u = Q * Binv_u;
  
  % Update inverse of the basis matrix
  Binv = Binv_u(:, 1:m);
end
endfunction
