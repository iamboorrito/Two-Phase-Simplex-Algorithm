# Two-Phase-Simplex-Algorithm
My implementation of the Two-Phase Simplex Algorithm for solving linear programs.

# Use

The main function is [x, z] = simplex(z0, c, A, b, ineqFlag) which solves the LP

        minimize z = c*x + z0
        subject to Ax=b.

where c is a row vector and the ineqFlag variable correspondes to -1 <-> Ax <= b, 0 <-> Ax = b, and likewise 1 <-> Ax >= b.

It does this by first adding artificial variables to form a basic solution to the
augmented problem, then attempts to drive them all to zero by elementary row operations
which exchange artificial variables for the original varialbes in the basis.

If there is a basic feasible solution, then the first phase will find it and then the second
phase prunes the next problem to remove redundant rows and finally calls the standard simplex
algorithm (rsimplex.m) to find the solution.

Ex. z0 = 0, c = [-2 -1 0 0], A = [2 1 1 0 ; 1 3 0 1], and b = [10 ; 15]

Enter [x, z] = simplex(0, c, A, b, 0)

where the last argument

Then the solution is given as x = [3 4 0 0]' with z = -10.

