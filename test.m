function x = test(A)
    % Multiply all rows of A by -1
    for i = 1:size(A,1)
        A(i,:) = -A(i,:);
    end
    % Return the first column of A
    x = A(:,1);
    return
endfunction