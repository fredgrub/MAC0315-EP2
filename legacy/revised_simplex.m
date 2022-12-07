function [z, x, pi, indices, exitflag] = RSM(A, b, c, m, n)
    % Solves min cx s.t. Ax=b, x>=0
    % exitflag is 0 if solved successfully, 1 if infeasible,and -1 if unbounded
    % Performs a Phase I procedure starting with an all artificial basis
    % and then calls function simplex
    
    % need to find Bmatrix, indices
    
    % set up the artificial indices
    indices = (n + 1):(n + m);
    % identity for the artificial
    IBmatrix = eye(m);
    % set up costs vector
    tempc = [zeros(n, 1); ones(m, 1)];
    cb = tempc(indices);
    
    [z, x, pi, indices, ~, IBmatrix] = simplex(A, b, tempc, m, n, IBmatrix, indices, 1);
    
    % in phase 2 we price all artificials as 0
    c = [c; zeros(m, 1)];

    % if a feasable solution has been found to phase 1
    if z == 0
        % we are left with the phase 1 complete, move to phase 2
        % in phase 2 we need to check for if each variable is artificial
        [z, x, pi, indices, exitflag, ~] = simplex(A, b, c, m, n, IBmatrix, indices, 2);
    else
        % problem is infeasable
        exitflag = 1;
    end
    
    % convert to column vector as required
    indices = indices';
end

function [z, x, pi, indices, exitflag, IBmatrix] = simplex(A, b, c, m, n, IBmatrix, indices, phase)
    % Solves min cx s.t. Ax=b, x>=0
    % starting with basic variables listed in vector indices
    % and basis matrix Bmatrix
    % returns optimal vector x and its value z, along with pi, and indices of basic variables
    
    exitflag = 0;
    xb = 0;
    cb = c(indices);
    
    while true
        pi = (cb' * IBmatrix)';
        
        [as, cs, s] = findenter(A, pi, c(1:n), indices);
        
        % check if the problem is optimal
        if s == 0
            break
        end
        
        % work out xb
        xb = IBmatrix * b;
        
        % leave returns the leaving column of the basis matrix
        % we pass in if each variable is artificial as well
        leave = findleave(IBmatrix, as, xb, phase, (indices > n)');

        % check if leave = 0 (unbounded)
        if leave == 0
            exitflag = -1;
            break;
        end
        
        [IBmatrix, indices, cb] = updateGJ(IBmatrix, indices, cb, cs, as, s, leave);
    end
    
    % calculate the final xb vector
    xb = IBmatrix * b;

    % get the new full x vector
    x = zeros(length(c), 1);
    x(indices) = xb;
    
    % calculate the new optimal cost
    z = c' * x;
    
    % shorten the x vector to disclude all artificals
    % indices could still contain artificals, but for the output of
    % the x vector, the artificals will always be value 0 so shouldn't 
    % need to be displayed
    x = x(1:n);
end

function [as, cs, s] = findenter(Amatrix, pi, c, indices)
    % Given the complete m by n matrix Amatrix,
    % the complete cost vector c with n components
    % the vector pi with m components
    % findenter finds the index of the entering variable and its column
    % It returns the column as, its cost coefficient cs, and index s
    % Returns s=0 if no entering variable can be found (i.e. optimal)
    
    % For the reduced costs, a reduced cost of 0 indicates that it is
    % already in the basis
    
    % an artificial variable is never non-basic
    non_basic = setdiff(1:length(c), indices);
    
    % if the if statement is never triggered, return that the problem is
    % optimal
    s = 0;
    as = [];
    cs = [];
    
    % pricing the columns individually
    for index = non_basic
        as = Amatrix(:, index);
        cs = c(index);
        rs = cs - (pi' * as)';
        
        if rs < -1e-8
            s = index;
            break;
        end
        
    end
end

function [leave] = findleave(IBmatrix, as, xb, phase, artificial)
    % Given entering column as and vector xb of basic variables
    % findleave finds a leaving column of basis matrix Bmatrix
    % It returns leave=0 if no column can be found (i.e.  unbounded)
    % artificial is logical array indicating if any variable is artificial
    
    % Make the 'important' vector
    important = IBmatrix * as;
    
    indices = 1:length(xb);
    
    % is any of the variables are artificial, we need to remove them first
    if phase == 2
        % get all of the artificial variables
        
        for i = 1:length(artificial)
            if artificial(i) && important(i) ~= 0
                leave = i;
                return;
            end
        end
        
        active = and((important > 0), ~artificial);
    else
        % Select all of the active columns which are valid for the ratio test
        active = (important > 0);
    end
    
    % Find the index of the leaving variable through the ratio test
    [~, leave_index] = min(xb(active) ./ important(active));
    
    % Check that there is a leaving variable
    if isempty(leave_index)
        % Problem is unbounded
        leave = 0;
    else
        % Get the correct leaving index from the active indices
        active_indices = indices(active);
        leave = active_indices(leave_index);
    end
end

function [IBmatrix, indices, cb] = updateGJ(IBmatrix, indices, cb, cs, as, s, leave)
    % In this function do the same as update, but update the inverse
    % of the B matrix using a gauss-jordan pivot.
    
    % costs basic
    % as is entering column
    % s is entering index
    % leave is leaving index
    
    important = IBmatrix * as;
    
    % we now need to pivot using the row leave
    stay = 1:length(indices) ~= leave;
    
    IBmatrix(stay, :) = IBmatrix(stay, :) - (important(stay) / important(leave)) * IBmatrix(leave, :);
    IBmatrix(leave, :) = IBmatrix(leave, :) / important(leave);
    
    indices(leave) = s;
    cb(leave) = cs;
end