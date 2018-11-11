#----------------------------------------
# File: my_model.mod (Wyndor_Glass_Co_linear_constraints)
# NEOS Server:  http://www.neos-server.org/neos/solvers
#----------------------------------------
reset;
param n; # Number of variables
param m; # Number of inequality constraints
param c{1..n}; # Objective function coefficients
param d{1..n}; # Objective function coefficients
param a{1..m,1..n}; # Inequality constraint coefficients
param b{1..m}; # Right-hand side inequality constraints
# Use either one or the other
data mymodel.dat;
#data mymodel_python_like.dat;
var x{1..n} >=0; #tons of water
var y{1..n} binary; # site binary



minimize z: # Objective function
sum {j in 1..n} (c[j]*x[j] +d[j]*y[j]);
subject to inequality {i in 1..m}: # Inequality constraints
sum {j in 1..n} a[i,j]*x[j] >= b[i];
subject to inequality2 {j in 1..n}: # Inequality constraints
x[j] <= 10000000000000*y[j];

data run.txt;
