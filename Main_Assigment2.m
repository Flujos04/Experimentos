clc;clear;close all

%% GRID PARAMETERS AND CONSTRUCTION %%

L=1;
n_cellsbump=10;
Nx=3*n_cellsbump+1;
Ny=60;

[x,y,y0,dx,dy]=mesher(L,Nx,Ny);

[x_center,y_center]=cell_center(Nx,Ny,x(1,:),y); %cell center

plot(x_center,y_center,'.')
hold off;

[A]=cell_area(Nx,Ny,x(1,:),y); %cell area (2D)

[Nx_N,Ny_N,Nx_E,Ny_E,Nx_S,Ny_S,Nx_W,Ny_W]=face_normals(Nx,Ny,x(1,:),y); %face normals

%% ENVIRONMENTAL PARAMETERS %%

p_ups=101300;          %Upstream pressure in Pa
T_ups=288;             %Upstream temperature in K
v_ups=0;               %Upstream vertical velocity
gamma=1.4;
R=287;
M_inf=0.1;           %Upstream Mach
a=sqrt(gamma*T_ups*R); %Upstream speed of sound
u_ups=a*M_inf; %Upstream horizontal velocity
%u_w=u_ups;
%v_w=0;
cv=R/(gamma-1);
cp=R*gamma/(gamma-1);

T_ups_s=T_ups*(1 + ((gamma-1)/2)*M_inf^2);                   %Upstream stagnation T in K
p_ups_s=p_ups*(1 + ((gamma-1)/2)*M_inf^2)^(gamma/(gamma-1)); %Upstream stagnation p in Pa
rho_ups=p_ups/(R*T_ups);                                       %Upstream density in kg/m^3
rho_ups_s=rho_ups*(p_ups/p_ups_s)^(gamma);                       %Upstream stagnation density in kg/m^3
E_ups=R/(gamma-1)*T_ups + (u_ups^2 + v_ups^2)/2;
H_ups=E_ups + p_ups/rho_ups;

%Prealocating variables
rho=zeros(Ny-1,Nx-1);
p=zeros(Ny-1,Nx-1);
T=zeros(Ny-1,Nx-1);
u=zeros(Ny-1,Nx-1);
v=zeros(Ny-1,Nx-1);
E=zeros(Ny-1,Nx-1);


%Initial conditions

rho(:,:)=rho_ups;
p(:,:)=p_ups;
T(:,:)=T_ups;
u(:,:)=u_ups;
v(:,:)=v_ups;
E(:,:)=R/(gamma-1)*T + (u.^2 + v.^2)/2;
H(:,:)=E + p./rho;

U(:,:,1) = rho;              
U(:,:,2) = rho.*u(:,:); 
U(:,:,3) = rho.*v(:,:);    
U(:,:,4) = rho.*E;
ntstep=10;

[U]=RungeKutta(U,ntstep,x,y,H,a,M_inf,u,rho,p,gamma,v,dx,dy,Nx,Ny,A,R,T,u_ups,v_ups,p_ups,Nx_E,Ny_E,Nx_W,Ny_W,Nx_S,Ny_S,Nx_N,Ny_N,cp,cv,T_ups_s,p_ups_s);
