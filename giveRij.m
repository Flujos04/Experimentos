function [Rij,FN,FE,FS,FW,f,g,dt]=giveRij(U,x,y,H,p,T,rho,u,v,p_ups,gamma,R,M_inf,a,u_ups,v_ups,Nx_E,Ny_E,Nx_W,Nx_S,Ny_W,Ny_S,Nx_N,Ny_N,Ny,Nx,CFL,k2,k4,cp,cv,T_ups_s,p_ups_s,dx,dy,rho_ups)

%%%%% WEST CONDITIONS (Inlet conditions) %%%%%
for j = 1:Ny-1        
    u_w(j,1) = (u_ups+u(j,1)+2*a/(gamma-1)-2*sqrt(gamma*p(j,1)/rho(j,1))/(gamma-1))/2;
    c_w(j,1) = (gamma-1)/2*(u_ups-u_w(j,1)+2*a/(gamma-1));
    T_w(j,1) = c_w(j,1)^2/(gamma*R);
    p_w(j,1) = (p_ups*(1/(R*T_w(j,1)*rho_ups))^gamma)^(1/(1-gamma));
    rho_w(j,1) = p_w(j,1)/(R*T_w(j,1));
    E_w(j,1) = p_w(j,1)/(gamma-1)/rho_w(j,1)+0.5*u_w(j,1)^2;
    H_w(j,1) = E_w(j,1)+p_w(j,1)/rho_w(j,1);    
    v_w(j,1) = v_ups;

    
    
    
    %     u(j,1) = (u_ups+u(j,1)+2*a/(gamma-1)-2*sqrt(gamma*p(j,1)/rho(j,1))/(gamma-1))/2;
%     c(j,1) = (gamma-1)/2*(u_ups-u(j,1)+2*a/(gamma-1));
%     T(j,1) = c(j,1)^2/(gamma*R);
%     p(j,1) = (p_ups*(1/(R*T(j,1)*rho_ups))^gamma)^(1/(1-gamma));
%     rho(j,1) = p(j,1)/(R*T(j,1));
%     E(j,1) = p(j,1)/(gamma-1)/rho(j,1)+0.5*u(j,1)^2;
%     H(j,1) = E(j,1)+p(j,1)/rho(j,1);    
%     v(j,1) = v_ups;

    
%     U(j,1,1)=rho_w(j,1);
%     U(j,1,2)=rho_w(j,1)*u_w(j,1);
%     U(j,1,3)=rho_w(j,1)*v_w(j,1);
%     U(j,1,4)=rho_w(j,1)*E_w(j,1);
       
end


% % T_w=T_ups_s-(gamma-1)./2.*(u_w.^2)./(gamma*R);
% % M_w=u_w./sqrt(gamma*R*T_w);
% % p_w=p_ups_s./(1+(gamma-1)./2*(M_w.^2))^(gamma/(gamma-1));
% % rho_w=p_w/(R*T_w); 
% % E_w=cv*T_w+(u_w.^2 + v_w.^2)/2;
% H_w=E_w+p_w./rho_w;

%%%%% EAST CONDITIONS (Outlet conditions)%%%%%
for j = 1:Ny-1
    p_e(j,1)=p_ups;
    rho_e(j,1) = (p_e(j,1)/p(j,Nx-1)*rho(j,Nx-1)^gamma)^(1/gamma);    
    T_e(j,1) = p_e(j,1)/(rho_e(j,1)*R);
    c_e(j,1) = sqrt(gamma*p_e(j,1)/rho_e(j,1));
    u_e(j,1) = u(j,Nx-1)+2*sqrt(gamma*p(j,Nx-1)/rho(j,Nx-1))/(gamma-1)-2*c_e(j,1)/(gamma-1);
    E_e(j,1) = p_e(j,1)/(gamma-1)/rho_e(j,1)+0.5*(u_e(j,1)^2+v(j,Nx-1)^2);
    H_e(j,1) = E_e(j,1)+p_e(j,1)/rho_e(j,1);
    
    v_e(j,1) = U(j,end,3)./U(j,end,1); 
    
%     U(j,Nx-1,1)=rho_e(j,1);
%     U(j,Nx-1,2)=rho_e(j,1)*u_e(j,1);
%     U(j,Nx-1,3)=rho_e(j,1)*v_e(j,1);
%     U(j,Nx-1,4)=rho_e(j,1)*E_e(j,1);
        
end

% p_e=p_ups;
% u_e = (U(:,end,2)./U(:,end,1));
% v_e = U(:,end,3)./U(:,end,1); 
% E_e = U(:,end,4)./U(:,end,1);
% H_e = E_e.*cp./cv;
% rho_e=U(:,end,1);


%%%%% NORTH CONDITIONS (Top wall)%%%%%
for j = 1:Nx-1    
    v_n(j,1)=0;

%     U(end,i,3)=rho(end,i)*v_n(i,1);
end

%%%%% SOUTH CONDITIONS (Bottom wall)%%%%%
for j = 1:Nx-1    
    v_s(j,1)=-u(1,j)*Nx_S(1,j)/Ny_S(1,j);

 %   U(1,i,3)=rho(end,i)*v(1,i);
end


% U(:,:,1)=rho;
% U(:,:,2)=rho.*u;
% U(:,:,3)=rho.*v;
% U(:,:,4)=rho.*E(j,1);




[f,g]=centre_fluxes_boundaries(Nx,Ny,rho,u,v,H,p,rho_w,u_w,p_w,v_w,H_w,rho_e,u_e,p_e,v_e,H_e,v_n,v_s);
[D,D_north,D_west,D_east,D_south] = artificial_dissipation (U,k4, k2, u, v, gamma, p, rho, Nx, Ny, x, y);
[FN, FE, FS, FW] = face_fluxes(f,g,Nx_N,Ny_N,Nx_E,Ny_E,Nx_S,Ny_S,Nx_W,Ny_W,Ny,Nx,D_north,D_west,D_east,D_south);
Rij=FN+FE+FS+FW;



u1=max(abs(u(:,:)+a));
u2=max(abs(u(:,:)-a));
v1=max(abs(v(:,:)+a));
v2=max(abs(v(:,:)-a));
umax=max([u1,u2]);
vmax=max([v1,v2]);
    
dt=CFL/((umax/dx) +(vmax/max(dy)));



end
