clc;clear;close all

L=1;
Nx=100;
Ny=100;

[x,y,y0]=mesher(L,Nx,Ny);

[x_center,y_center]=cell_center(Nx,Ny,x,y);

[A]=cell_area(Nx,Ny,x,y);
