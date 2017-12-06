function [ x, y, w ] = matrixGameSolver( A )

[m,n] = size(A);

inf_norm = norm(A, inf);

A = A + inf_norm;

c = ones(1,m);
b = ones(n, 1);

[x, w] = simplex(0, c, A', b, 1);
[y, ~] = simplex(0, -c, A, b, -1);

w = 1/w;
x = w*x;
y = w*y;
w = w - inf_norm;

end

