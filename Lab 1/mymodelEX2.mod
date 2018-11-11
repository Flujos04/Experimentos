#----------------------------------------
# File: my_model.mod (Wyndor_Glass_Co_linear_constraints)
# NEOS Server:  http://www.neos-server.org/neos/solvers
#----------------------------------------
reset;
param n; # Number of variables
param m; # Number of inequality constraints
param c{1..n}; # Objective function coefficients

param a{1..m,1..n}; # Inequality constraint coefficients

param b{1..m}; # Right-hand side inequality constraints

data mymodelEX2.dat;

var x{1..n} binary; # Variable definition
maximize z: # Objective function
sum {j in 1..n} c[j]*x[j];
subject to inequality {i in 1..m}: # Inequality constraints
sum {j in 1..n} a[i,j]*x[j] <= b[i];
data runEX2.txt;
