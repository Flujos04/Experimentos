reset;
param n;
param m;

param c{1..n};
param a{1..m,1..n};
param b{1..m};

#Data 
data ex3.dat;

#Variable definition

var x{1..n} >=4/12;

maximize z:
sum {j in 1..n} c[j]*x[j];

subject to inequality {i in 1..m}:
sum {j in 1..n} a[i,j]*x[j] <= b[i];

data runex3.txt;