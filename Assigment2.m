clc;clear;close all

L=1;
Nx=100;
Ny=100;
dx=3*L/(Nx-1);
x=0:dx:3*L;
R=(5/4 + 0.05)*L;
y0=zeros(1,Nx);
y0(1:Nx/3)=0; %BOTTOM LINE HEIGHT EQUATION 
y0((Nx/3 + dx):(Nx*2/3 - dx))=sqrt(R^2 - (x((Nx/3 + dx):(Nx*2/3 - dx))-1.5).^2) - (R - 0.1*L);%BOTTOM LINE HEIGHT EQUATION
y0(Nx*3*0.5:Nx)=0;%BOTTOM LINE HEIGHT EQUATION
dy=(L-y0)/(Ny-1);
y=0:dy:L;