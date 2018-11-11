reset;
param n; # Number of variables
param m; # Number of inequality constraints
param c{1..n}; # Objective function coefficients
param a{1..m,1..n}; # Inequality constraint coefficients
param b{1..m}; # Right-hand side inequality constraints
# Use either one or the other
data california.dat;
#data mymodel_python_like.dat;
var x{1..n} binary; # Variable definition
maximize z: # Objective function
sum {j in 1..n} c[j]*x[j];
subject to inequality {i in 1..m}: # Inequality constraints
sum {j in 1..n} a[i,j]*x[j] <= b[i];
data run.txt;
