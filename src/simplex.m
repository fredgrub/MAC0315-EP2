% Authors: Lucas de Sousa Rosa & Gabriel Schwartz
% Date: 05/12/2022
%
% Two phase revised simplex method for solving linear programming problems of
% the form
% 
%                           minimize c'x 
%                         subject to Ax = b
%                                    x >= 0
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

function [ind x d] = simplex(A,b,c,m,n)
    % phase I - find a feasible solution or the problem is infeasible

    % set up the artifical indices
    bind = (n + 1):(n + m);
    % B inverse
    Binv = eye(m);
    % set up cost vector
    tempc = [zeros(n,1); ones(m,1)];
    cB = tempc(bind);
    
    [z, x, ~, p, bind, ~, Binv] = RSM(A, b, tempc, m, n, Binv, bind, 1);

    % phase II - find the optimal solution or the problem is unbounded
    
    d = [];

    % change the cost vector making the artificial variables 0
    c = [c; zeros(m, 1)];

    % if a feasable solution has been found to phase 1,
    % we are left with the phase 1 complete, move to phase 2.
    if z == 0
        [z, x, d, p, bind, ind, ~] = RSM(A, b, c, m, n, Binv, bind, 2);
    else
        % problem is infeasable
        ind = 1;
    end
endfunction

function [z, x, d, p, bind, ind, Binv] = RSM(A, b, c, m, n, Binv, bind, phase)

    printf('**********\n');
    printf('* FASE %d *\n', phase);
    printf('**********\n\n');

    ind = 0;
    xb = 0;
    cb = c(bind);
    iter = 0;

    while true
        printf('==> Iteração %d\n\n', iter);
        iter = iter + 1;

        % calculate basic variables
        xb = Binv * b;

        % get the full x vector
        x = zeros(length(c), 1);
        x(bind) = xb;
        
        % calculate the optimal cost
        z = c' * x;

        printf('Variáveis básicas:\n');
        for i = 1:m
            printf('x%d = %f\n',bind(i), xb(i));
        endfor
        printf('\n');
        
        % Compute vector of simplex multipliers (p)
        p = (cb' * Binv)';

        printf('Custo atual da função objetivo: %f\n\n', z);

        [aj, rc, j] = findenter(A, p, c(1:n), bind);

        % check if the problem is optimal
        if j == 0
            break
        endif

        printf('Entra na base: %d\n\n', j);

        [l, d] = findleave(Binv, aj, xb, phase, (bind > n)', n, j, bind);

        % check if leave = 0 (unbounded)
        if l == 0
            ind = -1;
            break;
        endif

        printf('Sai da base: x%d\n\n',bind(l));

        % update inverse of the basis matrix using elementary row operations
        [Binv, bind, cb] = updateGJ(Binv, bind, cb, rc, aj, j, l);

    endwhile

    % calculate the final xb vector
    xb = Binv * b;

    % get the new full x vector
    x = zeros(length(c), 1);
    x(bind) = xb;
    
    % calculate the new optimal cost
    z = c' * x;
    
    % shorten the x vector to only the non-artifical variables
    x = x(1:n);

endfunction

function [aj, rc, j] = findenter(A, p, c, bind)
    % Find the entering variable using bland's rule
    %
    %   Input parameters:
    %
    % A - (m,n) constraint matrix 
    % p - (m,1) simplex multipliers vector
    % c - (n,1) cost vector
    % bind - (m,1) indices of current basic variables
    %
    %   Output parameters:
    % rc - is the reduced cost of the entering variable
    % aj - is the jth column (pivot column)
    % j - is the value of the entering variable
    %   j = 0 the problem is optimal (phase 1)

    % For the reduced costs, a reduced cost of 0 indicates that it is
    % already in the basis
    
    % an artificial variable is never non-basic
    non_basic = setdiff(1:length(c), bind);

    % if the if statement is never triggered, return that the problem is
    % optimal
    j = 0;
    aj = [];
    rc = [];
    
    % select entering variable using Bland's rule
    printf('Custos reduzidos (não-básica):\n');
    for index = non_basic
        aj = A(:, index);
        rc = c(index);
        rc = rc - (p' * aj)';

        printf('rc(x%d) = %f\n',index, rc);
        
        if rc < -1e-8
            j = index;
            break;
        endif
    endfor
    printf('\n');
endfunction

function [l, d] = findleave(Binv, aj, xb, phase, artificial, n, j, bind)
    % Find the leaving variable using bland's rule
    %
    %   Input parameters:
    %
    % Binv - (m,m) inverse of the basis matrix
    % aj - (m,1) entering column (pivot)
    % xb - (m,1) basic variables vector
    % phase - phase of the simplex method
    % artificial - (m,1) logical array indicating if any variable is artificial
    % n - number of variables (used for phase 2)
    % j - value of the entering variable (used for phase 2)
    %
    %   Output parameters:
    % l - is the value of the leaving variable
    %   l = 0 the problem is unbounded
    % d - (m,1) vector specifying the direction
    %   d is the direction on which the objective function goes to -infinity (when
    %   unbounded)

    d = [];
    
    % compute u = -db
    u = Binv * aj;

    indices = 1:length(xb);

    printf('Vetor de direção (apenas índices básicos):\n');
    for i = indices
        printf('u(x%d) = %f\n', bind(i), u(i));
    endfor
    printf('\n');

    % is any of the variables are artificial, we need to remove them first
    if phase == 2
        % get all of the artificial variables
        for i = 1:length(artificial)
            if artificial(i) && u(i) ~= 0
                l = i;
                return;
            endif
        endfor
        
        active = and((u > 0), ~artificial);
    else
        % Select all of the active columns which are valid for the ratio test
        active = (u > 0);
    endif

    % Find the index of the leaving variable through the ratio test
    [theta_star, leave_index] = min(xb(active) ./ u(active));
    
    printf('Theta*: %f\n\n', theta_star);

    % Check that there is a leaving variable
    if isempty(leave_index)
        % Problem is unbounded
        l = 0;
        d = zeros(1,n);
        d(j) = 1;
        d(indices) = -u;
    else
        % Get the correct leaving index from the active indices
        active_indices = indices(active);
        l = active_indices(leave_index);
    endif
endfunction

function [Binv, bind, cb] = updateGJ(Binv, bind, cb, rc, aj, j, l)
    % update the inverse of the B matrix using elementary row operations.
    %
    %  Input parameters:
    % Binv - (m,m) inverse of the basis matrix
    % bind - (m,1) indices of current basic variables
    % cb - (m,1) cost vector of the basic variables
    % rc - is the reduced cost of the entering variable
    % aj - is the jth column (pivot column)
    % j - is the value of the entering variable
    % l - is the value of the leaving variable
    %
    %  Output parameters:
    % Binv - (m,m) inverse of the basis matrix after the update
    % bind - (m,1) indices of current basic variables after the update
    % cb - (m,1) cost vector of the basic variables after the update

    u = Binv * aj;

    % pivot using the row leave
    stay = 1:length(bind) ~= l;

    % update the inverse of the basis matrix
    Binv(stay, :) = Binv(stay, :) - (u(stay) / u(l)) * Binv(l, :);
    Binv(l, :) = Binv(l, :) / u(l);

    bind(l) = j;
    cb(l) = rc;
endfunction