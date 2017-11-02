


val = 10000

m = 10;
n = 15;

while val ~= 1
    
    %A = rand(m, n) ;
    b = rand(m,1);
    c = rand(1, n);
    [~, ~, val] = linprog(c, A, b);
    
end

c

A

b

