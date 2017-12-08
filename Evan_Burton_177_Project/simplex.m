function [ x_opt, z_opt ] = simplex(z0, c, A, b, ineqFlag)
% simplex(z0, c, A, b, ineqFlag) solves the linear program given by
%
%        minimize z = c*x + z0
%        subject to Ax (<= / = / >=) b.
%       ineqFlag:=      -1   0    1
%
% where c is a row vector.
%
%-----------------------------------------------
% Example:  c = [-2 -3 -3]
%           A = [3 2 0 ; -1 1 4 ; 2 -2 5]           
%           b = [60 ; 10 ; 50]
%   
% For minimize c*x + 10 subject to Ax <= b:
% 
% [x, z] = simplex(10, c, A, b, -1)
%  x = [8, 18, 0]' and z = -60
%
% For minimize c*x + 10 subject to Ax = b:
%
% [x, z] = simplex(10, c, A, b, 0)
%------------------------------------------------
% Example: c = [-1 -2 -3]
%          A = [1 1 -1 ; -2 1 2 ; 1 -1 0 ; 0 1 1]
%          b = [1 ; 5 ; 4 ; 5]
% 
% [x, z] = simplex(0, c, A, b, -1)
%  x = [4, 0, 5]' and z = -19
%------------------------------------------------
% 
% For the included redundant system:
% A = csvread('redA.csv')
% b = csvread('redb.csv')
% c = csvread('redc.csv')
%
% [x, z] = simplex(0, c, A, b, 0)
% and x4 = 10, xk = 0 for k ~= 4.

% This is for record keeping purposes, need -z0 for tableau
z0 = -z0;
[m, n] = size(A);

if ineqFlag == -1
    [ x_opt, z_opt, ~, ~] = rsimplex(z0,[c, zeros(1, m)], [A, eye(m)], b, n+1:n+m);
    x_opt = x_opt(1:n);
    return;
else
if ineqFlag == 1
    A = [A, -eye(m)];
    %w = [-sum(A, 1) zeros(1,m)];
    c = [c, zeros(1, m)];
    n = n+m;
else
    %w = [-sum(A, 1), zeros(1,m)];
end

% New objective
w = [-sum(A, 1), zeros(1,m)];

%[[A eye(m)], b; w, -sum(b)]

% Solve augmented problem
[ x_opt, z_opt, CARRY, basis ] = rsimplex(-sum(b), w, [A, eye(m)], b, []);

%ABS_ZOPT = abs(z_opt)
if abs(z_opt) > 10^(-10)
    x_opt = NaN;
    z_opt = NaN;
    MSG = sprintf('%s', 'Infeasible solution')
    return
end

%BFS = x_opt
%MSG = sprintf('%s', 'FEASIBLE SOLN FOUND, BEGIN PHASE TWO')

A_star = CARRY(2:m+1, 2:m+1)*A;

% Remove redundant rows, found where basis(i) > n
% since artificial variable will have index > n
rows_to_remove = sum(basis > n);
rows_to_keep = m-rows_to_remove;
rows = zeros(1, rows_to_keep-1);

if rows_to_remove > 0
    k = 1;
    for i = 1:m
        
        if basis(i) < n && k <= rows_to_keep
            
            rows(k) = i;
            
            k = k+1;
        end
        
    end
    
    A_star = A_star(rows, :);
    
    basis = basis(rows);
    CARRY = CARRY([1, 1+rows], 1:m);
    
    m = size(A_star, 1);
end

% Calculate the new entries of the table
c_star = c-c(basis)*A_star;
b_star = CARRY(2:m+1, 1);
z0_star = z0 - c(basis)*b_star;

%[ A_star b_star; c_star z0_star]

%pause

[x_opt, z_opt, ~] = rsimplex(z0_star, c_star, A_star, b_star, basis);

% Only give values for original variables
if ineqFlag == 1 && abs(z_opt) < inf
    x_opt = x_opt(1:n-m);
end

end

