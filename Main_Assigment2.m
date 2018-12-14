clc;clear;close all

%% GRID PARAMETERS AND CONSTRUCTION %%

L=1;
n_cellsbump=10;
Nx=3*n_cellsbump+1;
Ny=60;

[x,y,y0]=mesher(L,Nx,Ny); %mesh

[x_center,y_center]=cell_center(Nx,Ny,x,y); %cell center

plot(x_center,y_center,'.')
hold off;

[A]=cell_area(Nx,Ny,x,y); %cell area (2D)

[Nx_N,Ny_N,Nx_E,Ny_E,Nx_S,Ny_S,Nx_W,Ny_W]=face_normals(Nx,Ny,x,y); %face normals

%% ENVIRONMENTAL PARAMETERS %%

p_ups=101300;          %Upstream pressure in Pa
T_ups=288;             %Upstream temperature in K
v_ups=0;               %Upstream vertical velocity
gamma=1.4;
R=287;
M_inf=0.1;           %Upstream Mach
a=sqrt(gamma*T_w*R); %Upstream speed of sound
u_ups=a*M_inf; %Upstream horizontal velocity
u_w=u_ups;
v_w=0;
cp=R./(gamma-1);
cv=R*gamma./(gamma-1);

T_w_s=T_ups*(1 + ((gamma-1)/2)*M_inf^2);                   %Upstream stagnation T in K
p_w_s=p_ups*(1 + ((gamma-1)/2)*M_inf^2)^(gamma/(gamma-1)); %Upstream stagnation p in Pa
rho_ups=p_ups/(R*T_ups);                                       %Upstream density in kg/m^3
rho_w_s=rho_wups(p_w/p_w_s)^(gamma);                       %Upstream stagnation density in kg/m^3
E=R/(gamma-1)*T_ups + (u_w^2 + v_w^2)/2;
H=E + p_w/rho_ups;

%Prealocating variables
rho=zeros(Ny-1,Nx-1);
p=zeros(Ny-1,Nx-1);
T=zeros(Ny-1,Nx-1);
u=zeros(Ny-1,Nx-1);
v=zeros(Ny-1,Nx-1);
E=zeros(Ny-1,Nx-1);


%Initial conditions

rho(:,:)=rho_w;
p(:,:)=p_w;
T(:,:)=T_w;
u(:,:)=u_w;
v(:,:)=v_w;
E(:,:)=R/(gamma-1)*T + (u.^2 + v.^2)/2;

U(:,:,1) = rho;              
U(:,:,2) = rho.*u(:,:); 
U(:,:,3) = rho.*v(:,:);    
U(:,:,4) = rho.*E;
ntstep=10;
[U]=RungeKutta(U,ntstep,a,u,rho,p,gamma,v,dx,dy,Nx,Ny,A,R,T,Nx_E,Ny_E,Nx_W,Nx_S,Ny_w,Ny_S,Nx_N,Ny_N,T_w_s,p_w_s,rho_w_s)

