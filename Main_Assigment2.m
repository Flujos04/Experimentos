clc;clear;close all

L=1;
n_cellsbump=100;
Nx=3*n_cellsbump+1;
Ny=60;

[x,y,y0]=mesher(L,Nx,Ny); %mesh

[x_center,y_center]=cell_center(Nx,Ny,x,y); %cell center

plot(x_center,y_center,'.')
hold off;

[A]=cell_area(Nx,Ny,x,y); %cell area (2D)

[dS]=face_normals(Nx,Ny,x,y); %face normals
