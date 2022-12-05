function [ind x d] = simplex(A,b,c,m,n)
    % Two phase revised simplex method for solving linear programming problems
    % of the form
    % 
    %                         min        c'x 
    %                         subject to Ax = b, 
    %                                    x >= 0.
    %
    % Input: A, b, c, m, n
    % A = m x n constraint matrix 
    % b = m x 1 right-hand side vector
    % c = n x 1 cost vector
    % m = number of constraints
    % n = number of variables
    %
    % Output: ind, x, d
    % ind = -1 if the problem is unbounded
    % ind = 1 if the problem is infeasible
    % ind = 0 if the problem has optimal solution
    % x = n x 1 vector
    % d = m x 1 vector
    %   x is the last basic feasible solution found. If the problem is unbounded,
    %   d is the direction on which the objective function goes to -infinity.

    % Solve the auxiliary problem (phase 1)
    [ind,x,bind,Binv] = auxiliary_problem(A,b,c,m,n);

    if ind == 0
        % The auxiliary problem has optimal solution. Now we solve the original
        % problem (phase 2).
        [ind,x,~,~] = solver(A,b,c,m,n,x,bind,Binv,2);

        if ind == 0
            printf('Solução ótima encontrada. Valor da função objetivo: %f\n',c'*x);
            for i = 1:n
                printf('x%d = %f\n',i,x(i));
                return
            end
        elseif ind == -1
            printf('Solução Ilimitada.\n');
            printf('Direção ilimitada:\n');
            for i = 1:n
                printf('d%d = %f\n',i,x(i));
            end
            return
        end
    elseif ind == 1
        printf('The problem is infeasible');
        return
    end

endfunction

function [ind,x,bind,Binv] = auxiliary_problem(A,b,c,m,n)
    % Solve the auxiliary problem to determine an initial basic feasible solution.
    % The auxiliary problem is of the form
    %                 min        0'x + 1'y
    %                 subject to Ax = b 
    %                            x >= 0
    %                            y >= 0
    %
    % Output: ind, v
    % ind = 1 if the optimal cost in the auxiliary problem is positive (the original
    % problem is infeasible)
    % ind = 0 if the optimal cost in the auxiliary problem  is zero.
    %   In the latter case, v is the initial basic feasible solution.

    % Check if b >= 0 and if not, multiply the corresponding row of A and b by -1
    for i = 1:m
        if b(i) < 0
            A(i,:) = -A(i,:);
            b(i) = -b(i);
        end
    end

    % Add artificial variables to the problem and solve it
    A = [A eye(m)];
    c = [zeros(n,1); ones(m,1)];
    n = n + m;

    % Define x as the initial basic feasible solution
    x = zeros(n,1);
    x(n-m+1:n) = b;

    % Define bind as the set of indices of the basic variables
    bind = n-m+1:n;

    % Compute Binv, the inverse of the basis matrix B
    Binv = inv(A(:,bind));

    % Solve the auxiliary problem
    [ind,x,bind,Binv] = solver(A,b,c,m,n,x,bind,Binv,1)
    return
endfunction