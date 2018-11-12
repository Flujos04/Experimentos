#----------------------------------------
# File: my_model.mod (Wyndor_Glass_Co_linear_constraints)
# NEOS Server:  http://www.neos-server.org/neos/solvers
#----------------------------------------
reset;
param n; # Number of variables
param m; # Number of inequality constraints

param c{1..n}; # Objective function coefficients

param a{1..m,1..n}; # Inequality constraint coefficients
#param aeq{1..l,1..n}; # equality constraint coefficients (they do not apply for this problem)

param b{1..m}; # Right-hand side inequality constraints
#param beq{1..l}; # Right-hand side equality constraints (they do not apply for this problem)

#data mymodel.dat;
data mymodel5.dat;

var x{1..n} binary; # Variable definition

minimize z: # Objective function
sum {j in 1..n} c[j]*x[j];

subject to inequality {i in 1..m}: # Inequality constraints
sum {j in 1..n} a[i,j]*x[j] >= b[i];

subject to equality: # Inequality constraints
sum {j in 1..n} x[j] = 3;

#equality {k in 1..l}: # Equality constraints (they do not apply for this problem)
#sum {j in 1..n} aeq[k,j]*x[j] = beq[k]; (they do not apply for this problem)

data run5.txt;
