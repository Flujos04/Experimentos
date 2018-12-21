function [FN, FE, FS, FW] = face_fluxes(f,g,Nx_N,Ny_N,Nx_E,Ny_E,Nx_S,Ny_S,Nx_W,Ny_W,Ny,Nx,D_north,D_west,D_east,D_south,u_ups,a,rho_ups,p_ups,v_ups,u,rho,p,gamma,R,v)
         
F_ew = zeros(Ny-1,Nx,4); %fluxes in east and west faces
G_ew = zeros(Ny-1,Nx,4); %fluxes in east and west faces
F_sn = zeros(Ny,Nx-1,4); %fluxes in south and north faces
G_sn = zeros(Ny,Nx-1,4); %fluxes in south and north faces


%Calculation of fluxes in east and west faces
for j = 1:Ny-1
    for i = 1:Nx-2
        F_ew(j,i+1,:) = (f(j,i,:)+f(j,i+1,:))/2;
        G_ew(j,i+1,:) = (g(j,i,:)+g(j,i+1,:))/2;        
    end
end

%Calculation of fluxes in north and south faces
for j = 1:Ny-2
    for i = 1:Nx-1
        F_sn(j+1,i,:) = (f(j,i,:)+f(j+1,i,:))/2;
        G_sn(j+1,i,:) = (g(j,i,:)+g(j+1,i,:))/2;
    end
end


%% BOUNDARY CONDITIONS
% Inlet boundary conditions
    for j = 1:Ny-1
        u_inlet(j,1) = (u_ups+u(j,1)+2*a/(gamma-1)-2*sqrt(gamma*p(j,1)/rho(j,1))/(gamma-1))/2;
        c_inlet(j,1) = (gamma-1)/2*(u_ups-u_inlet(j,1)+2*a/(gamma-1));
        T_inlet(j,1) = c_inlet(j,1)^2/(gamma*R);
        p_inlet(j,1) = (p_ups*(1/(R*T_inlet(j,1)*rho_ups))^gamma)^(1/(1-gamma));
        rho_inlet(j,1) = p_inlet(j,1)/(R*T_inlet(j,1));
        E_inlet(j,1) = p_inlet(j,1)/(gamma-1)/rho_inlet(j,1)+0.5*u_inlet(j,1)^2;
        H_inlet(j,1) = E_inlet(j,1)+p_inlet(j,1)/rho_inlet(j,1);
        v_inlet(j,i) = v_ups;
    end

    for j = 1:Ny-1
        F_ew(j,1,1) = rho_inlet(j,1)*u_inlet(j,1);
        F_ew(j,1,2) = rho_inlet(j,1)*u_inlet(j,1)^2+p_inlet(j,1);
        F_ew(j,1,3) = rho_inlet(j,1)*u_inlet(j,1)*v_inlet(j,1); %=0
        F_ew(j,1,4) = rho_inlet(j,1)*u_inlet(j,1)*H_inlet(j,1);
        
        G_ew(j,1,1) = rho_inlet(j,1)*v_inlet(j,1);%=0
        G_ew(j,1,2) = rho_inlet(j,1)*u_inlet(j,1)*v_inlet(j,1);%=0
        G_ew(j,1,3) = rho_inlet(j,1)*v_inlet(j,1)^2+p_inlet(j,1);%=p_inlet
        G_ew(j,1,4) = rho_inlet(j,1)*v_inlet(j,1)*H_inlet(j,1); %=0
        %Then, G will be removed because it will be multiplied by Ny (which is 0)
    end

    
% Outlet boundary conditions
    for j = 1:Ny-1
        p_outlet(j,1)=p_ups;
        rho_outlet(j,1) = (p_outlet(j,1)/p(j,Nx-1)*rho(j,Nx-1)^gamma)^(1/gamma);
        T_outlet(j,1) = p_outlet(j,1)/(rho_outlet(j,1)*R);
        c_outlet(j,1) = sqrt(gamma*p_outlet(j,1)/rho_outlet(j,1));
        u_outlet(j,1) = u(j,Nx-1)+2*sqrt(gamma*p(j,Nx-1)/rho(j,Nx-1))/(gamma-1)-2*c_outlet(j,1)/(gamma-1);
        E_outlet(j,1) = p_outlet(j,1)/(gamma-1)/rho_outlet(j,1)+0.5*(u_outlet(j,1)^2+v(j,Nx-1)^2);
        H_outlet(j,1) = E_outlet(j,1)+p_outlet(j,1)/rho_outlet(j,1);
    end

    for j = 1:Ny-1
        F_ew(j,Nx,1) = rho_outlet(j,1)*u_outlet(j,1);
        F_ew(j,Nx,2) = rho_outlet(j,1)*u_outlet(j,1)^2+p_outlet(j,1);
        F_ew(j,Nx,3) = rho_outlet(j,1)*u_outlet(j,1)*v(j,Nx-1);
        F_ew(j,Nx,4) = rho_outlet(j,1)*u_outlet(j,1)*H_outlet(j,1);
        
        G_ew(j,Nx,1) = rho_outlet(j,1)*v(j,Nx-1);
        G_ew(j,Nx,2) = rho_outlet(j,1)*u_outlet(j,1)*v(j,Nx-1);
        G_ew(j,Nx,3) = rho_outlet(j,1)*v(j,Nx-1)^2+p_outlet(j,1);
        G_ew(j,Nx,4) = rho_outlet(j,1)*v(j,Nx-1)^2*H_outlet(j,1);    
        %Then, G will be removed because it will be multiplied by Ny (which is 0)
    end
    

    
% Lower wall boundary conditions
    for i = 1:Nx-1
        F_sn(1,i,2)=p(1,i);
        G_sn(1,i,3)=p(1,i);
    end

% Upper wall boundary conditions
    for i = 1:Nx-1
        G_sn(Ny,i,3) = p(Ny-1,i);
    end


FN = zeros((Ny-1),(Nx-1),4);
FE = zeros((Ny-1),(Nx-1),4);
FS = zeros((Ny-1),(Nx-1),4);
FW = zeros((Ny-1),(Nx-1),4);

for i=1:Nx-1
    for j=1:Ny-1
        FN(j,i,:) =F_sn(j+1,i,:)*Nx_N(j,i) + G_sn(j+1,i,:)*Ny_N(j,i)-D_north(j,i,:);
        FS(j,i,:) =F_sn(j,i,:)*Nx_S(j,i) + G_sn(j,i,:)*Ny_S(j,i)+D_south(j,i,:);
        FE(j,i,:) =F_ew(j,i+1,:)*Nx_E(j,i) + G_ew(j,i+1,:)*Ny_E(j,i)-D_east(j,i,:);
        FW(j,i,:) =F_ew(j,i,:)*Nx_W(j,i) + G_ew(j,i,:)*Ny_W(j,i)+D_west(j,i,:);
        
    end
end

end
