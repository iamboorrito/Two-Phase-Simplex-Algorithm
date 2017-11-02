function [ x_opt, z_opt ] = simplex(z0, c, A, b)
%Solves the linear program given by 
%
%        minimize z = c*x + z0
%        subject to Ax=b.
%
% where c is a row vector.
[m, n] = size(A);

% New objective
w = -sum(A, 1);

MSG = sprintf('%s', 'BEGIN PHASE ONE')
[ x_opt, z_opt, CARRY, basis ] = rsimplex(-sum(b), w, A, b, []);

if abs(z_opt) > 10^(-14)
    %b_opt = NaN;
    z_opt = NaN;
    return
end

BFS_TRANSPOSE = x_opt'
MSG = sprintf('%s', 'FEASIBLE SOLN FOUND, BEGIN PHASE TWO')

A_star = CARRY(2:m+1, 2:m+1)*A;

%Remove redundant rows, found where basis(i) > n
rows_to_remove = sum(basis > n);
rows_to_keep = m-rows_to_remove;
rows = zeros(1, rows_to_keep);

if rows_to_remove > 0
    k = 1;
    for i = 1:m

        if basis(i) < n && k <= rows_to_keep

            rows(k) = i;

            k = k+1;
        end

    end

    %rows
    A_star = A_star(rows, :);
    basis = basis(rows);
    CARRY = CARRY([1, 1+rows], 1:m);

    m = size(A_star, 1);
end

% Calculate the new entries of the table
c_star = c-c(basis)*A_star;
b_star = CARRY(2:m+1, 1);
z0_star = z0 - c(basis)*b_star;

%pause

[x_opt, z_opt, ~] = rsimplex(z0_star, c_star, A_star, b_star, basis);

end

