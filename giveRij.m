function [Rij,FN,FE,FS,FW,f,g,dt]=giveRij(U,x,y,H,p,T,rho,u,v,p_ups,gamma,R,M_inf,a,u_ups,v_ups,Nx_E,Ny_E,Nx_W,Nx_S,Ny_W,Ny_S,Nx_N,Ny_N,Ny,Nx,CFL,k2,k4,cp,cv,T_ups_s,p_ups_s,dx,dy)

%%%%% WEST CONDITIONS %%%%%


u_w=u_ups;
v_w=v_ups;

T_w=T_ups_s-(gamma-1)./2.*(u_w.^2)./(gamma*R);
M_w=u_w./sqrt(gamma*R*T_w);
p_w=p_ups_s./(1+(gamma-1)./2*(M_w.^2))^(gamma/(gamma-1));
rho_w=p_w/(R*T_w); 
E_w=cv*T_w+(u_w.^2 + v_w.^2)/2;
H_w=E_w+p_w./rho_w;

%%%%% EAST CONDITIONS %%%%%

p_e=p_ups;
u_e = (U(:,end,2)./U(:,end,1));
v_e = U(:,end,3)./U(:,end,1); 
E_e = U(:,end,4)./U(:,end,1);
H_e = E_e.*cp./cv;
rho_e=U(:,end,1);

% [f,g]=centre_fluxes(Nx,Ny,rho_w,p_w,H_w,rho_w,p_w,H_w,u_w,v_w,rho_e,p_e,H_e,u_e,v_e);
[f,g,dt]=centre_fluxes_boundaries(Nx,Ny,rho,u,v,H,p,rho_w,p_w,H_w,u_w,v_w,rho_e,p_e,H_e,u_e,v_e,a,CFL,dx,dy);
[D,D_north,D_west,D_east,D_south] = artificial_dissipation (U,k4, k2, u, v, gamma, p, rho, Nx, Ny, x, y);
[FN, FE, FS, FW] = face_fluxes(f,g,Nx_N,Ny_N,Nx_E,Ny_E,Nx_S,Ny_S,Nx_W,Ny_W,Ny,Nx,D_north,D_west,D_east,D_south);
Rij=FN+FE+FS+FW;

end

