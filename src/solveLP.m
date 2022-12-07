function solveLP()
    [A,b,c,m,n] = get_problem();
    [ind x d] = simplex(A,b,c,m,n);

    if (ind == -1)
        printf('### Problema ilimitado. Direção de custo -INFINITO:\n');
        for i = 1:n
            printf('x%d = %f\n',i,d(i));
        endfor
    elseif (ind == 0)
        printf('### Solução ótima encontrada. Valor da função objetivo: %f\n',c'*x);
        for i = 1:n
            printf('x%d = %f\n',i,x(i));
        endfor
    else
        printf('### Problema inviável.\n');
    endif
endfunction

function [A,b,c,m,n] = get_problem()
    % Unbounded (ex. 1 da P2)
    A = [2 1 -2 1 0 0; 4 1 2 0 -1 0; 2 3 -1 0 0 -1];
    b = [8; 2; 4];
    c = [-2; 1; -1; 0; 0; 0];
    m = 3;
    n = 6;

    % Feasible (ex. 3.17)
    % A = [1 3 0 4 1; 1 2 0 -3 1; -1 -4 3 0 0];
    % b = [2; 2; 1];
    % c = [2; 3; 3; 1; -2];
    % m = 3;
    % n = 5;

    % Infeasible (inventado)
    % A = [1 0 -1 0 0; 0 1 0 -1 0; 1 1 0 0 1];
    % b = [6; 6; 11];
    % c = [1; 1; 0; 0; 0];
    % m = 3;
    % n = 5;
endfunction