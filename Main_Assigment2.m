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

%% ENVIRONMENTAL PARAMETERS %%

p_w=101300;          %Upstream pressure in Pa
T_w=288;             %Upstream temperature in K
v_w=0;               %Upstream vertical velocity
gamma=1.4;
R=287;
M_inf=0.1;           %Upstream Mach
a=sqrt(gamma*T_w*R); %Upstream speed of sound
u_w=a*M_inf;         %Upstream horizontal velocity

T_w_s=T_w*(1 + ((gamma-1)/2)*M_inf^2);                   %Upstream stagnation T in K
p_w_s=p_w*(1 + ((gamma-1)/2)*M_inf^2)^(gamma/(gamma-1)); %Upstream stagnation p in Pa
rho_w=p_w/(R*T_w);                                       %Upstream density in kg/m^3
rho_w_s=rho_w*(p_w/p_w_s)^(gamma);                       %Upstream stagnation density in kg/m^3
E=R/(gamma-1)*T_w + (u_w^2 + v_w^2)/2;
H=E + p_w/rho_w;

%Prealocating variables (Not counting the boundaries)

rho=zeros(Ny-1,Nx-1);
p=zeros(Ny-1,Nx-1);
T=zeros(Ny-1,Nx-1);
u=zeros(Ny-1,Nx-1);
v=zeros(Ny-1,Nx-1);
E=zeros(Ny-1,Nx-1);


%Initial conditions (Not counting the boundaries)

rho0=rho*rho_w;
T0=T*T_w;
u0=u*u_w;
v0=v*v_q;
E0=R/(gamma-1)*T0 + (u0^2 + v0^2)/2;

