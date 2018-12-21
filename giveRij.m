function [Rij,FN,FE,FS,FW,f,g,dt]=giveRij(U,x,y,H,p,T,rho,u,v,p_ups,gamma,R,M_inf,a,u_ups,v_ups,Nx_E,Ny_E,Nx_W,Nx_S,Ny_W,Ny_S,Nx_N,Ny_N,Ny,Nx,CFL,k2,k4,cp,cv,T_ups_s,p_ups_s,dx,dy,rho_ups)


[f,g] = centre_fluxes_boundaries(Nx,Ny,rho,u,v,H,p);
[D,D_north,D_west,D_east,D_south] = artificial_dissipation (U,k4, k2, u, v, gamma, p, rho, Nx, Ny, x, y);
[FN, FE, FS, FW] = face_fluxes(f,g,Nx_N,Ny_N,Nx_E,Ny_E,Nx_S,Ny_S,Nx_W,Ny_W,Ny,Nx,D_north,D_west,D_east,D_south,u_ups,a,rho_ups,p_ups,v_ups,u,rho,p,gamma,R,v);
Rij=FN+FE+FS+FW;


u1=max(abs(u(:,:)+a));
u2=max(abs(u(:,:)-a));
v1=max(abs(v(:,:)+a));
v2=max(abs(v(:,:)-a));
umax=max([u1,u2]);
vmax=max([v1,v2]);
    
dt=CFL/((umax/dx) +(vmax/max(dy)));



end
