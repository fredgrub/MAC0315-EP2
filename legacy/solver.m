function [ind, v, bind, Binv] = solver(A,b,c,m,n,x,bind,Binv,phase)
  % Revised Simplex solver for solving the linear programming problem
  % of the form
  % 
  %                         min        c'x 
  %                         subject to Ax = b, 
  %                                    x >= 0.
  %
  % Input: A, b, c, m, n, x, bind, Binv, phase
  % A = m x n constraint matrix
  % b = m x 1 right-hand side vector
  % c = n x 1 cost vector
  % m = number of constraints
  % n = number of variables
  % x = n x 1 initial basic feasible solution
  % bind = m x 1 vector of indices of basic variables
  % Binv = m x m inverse of the basis matrix
  % phase = current phase of the algorithm (1 or 2)
  %
  % Output: ind, v, tableau
  % ind = 0 if the solution is optimal
  % ind = -1 if the solution is unbounded
  % v = optimal solution if ind = 0
  % v = unbounded direction if ind = -1

  printf('Fase %d\n', phase)
  iter = 0;

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
      % on phase 1, check if the optimal cost is positive (infeasible)
      % and the algorithm terminates. Otherwise, if the optimal cost is
      % zero and no artificial variables are in the basis, the algorithm
      % continues to phase 2. Otherwise, the artificial variables and the
      % corresponding constraints are removed from the problem and the 
      % algorithm continues to phase 2. If any basic variable is artificial,
      % the phase 1 algorithm continues.
      if phase == 1
        if c'*x > 0
          ind = 1;
          v = [];
          return;
        elseif c'*x == 0 && isempty(intersect(bind,n - m + 1:n))
          ind = 0;
          % remove artificial variables from basic solution
          x = x(1:n);
          v = x;
          return
        end
      else
        % Optimal solution found (phase 2)
        ind = 0;
        v = x;
        return;
      end
    end

    %%%%%%%%%%%%%%%%%%%%%% UNDER CONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%

    % Select entering variable using Bland's rule
    for i = 1:length(nonbind)
      if rc(i) < 0
        j = nonbind(i);
        break
      end
    end
  
    % Display entering variable
    printf('Entra na base: x%d\n\n',j);

    % Compute direction vector
    u = Binv*A(:,j);

    % Check unboundedness
    if phase == 2 && all(u <= 0)
        ind = -1;
        % Calculate unbounded direction
        d = zeros(n,1);
        d(j) = 1;
        d(bind) = -u;
        v = d;
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

    if iter > 3
      printf('Número máximo de iterações atingido.\n');
      return
    end
  end
endfunction