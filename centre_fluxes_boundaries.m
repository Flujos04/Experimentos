function [f,g]=centre_fluxes_boundaries(Nx,Ny,rho,ux,uy,H,p,rho_w,p_w,H_w,u_w,v_w,rho_e,p_e,H_e,u_e,v_e)

%We are going to work with F and G in matrixes. Then, they will be converted to vectors for simplicity
f = zeros((Ny-1),(Nx-1),4);
g = zeros((Ny-1),(Nx-1),4);


for i=2:(Nx-2)
    for j=2:(Ny-2)
        
        %%%%% INNNER CELLS %%%%%
        f(j,i,1)=rho(j,i)*ux(j,i);
        f(j,i,2)=rho(j,i)*ux(j,i)^2+p(j,i);
        f(j,i,3)=rho(j,i)*ux(j,i)*uy(j,i);
        f(j,i,4)=rho(j,i)*ux(j,i)*H(j,i);
        
        g(j,i,1)=rho(j,i)*uy(j,i);
        g(j,i,2)=rho(j,i)*ux(j,i)*uy(j,i);
        g(j,i,3)=rho(j,i)*uy(j,i)^2+p(j,i);
        g(j,i,4)=rho(j,i)*uy(j,i)*H(j,i);
        
        %%%%% WEST CELLS %%%%%
        f(j,1,1)=rho_w*u_w;
        f(j,1,2)=rho_w*u_w^2+p_w;
        f(j,1,3)=rho_w*u_w*v_w;
        f(j,1,4)=rho_w*u_w*H_w;
        
        g(j,1,1)=rho_w*v_w;
        g(j,1,2)=rho_w*u_w*v_w;
        g(j,1,3)=rho_w*v_w^2+p_w;
        g(j,1,4)=rho_w*v_w*H_w;
        
        %%%%% EAST CELLS %%%%%
        f(j,1,1)=rho_e(j)*u_e(j);
        f(j,1,2)=rho_e(j)*u_e(j)^2+p_e;
        f(j,1,3)=rho_e(j)*u_e(j)*v_e(j);
        f(j,1,4)=rho_e(j)*u_e(j)*H_e(j);
        
        g(j,1,1)=rho_e(j)*v_e(j);
        g(j,1,2)=rho_e(j)*u_e(j)*v_e(j);
        g(j,1,3)=rho_e(j)*v_e(j)^2+p_e;
        g(j,1,4)=rho_e(j)*v_e(j)*H_e(j);
        
    end
    
        %%%%% NORTH CELLS %%%%%
        f(end,i,1)=rho(end,i)*ux(end,i);
        f(end,i,2)=rho(end,i)*ux(end,i)^2+p(end,i);
        f(end,i,3)=rho(end,i)*ux(end,i)*uy(end,i);
        f(end,i,4)=rho(end,i)*ux(end,i)*H(end,i);
        
        g(end,i,1)=rho(end,i)*uy(end,i);
        g(end,i,2)=rho(end,i)*ux(end,i)*uy(end,i);
        g(end,i,3)=rho(end,i)*uy(end,i)^2+p(end,i);
        g(end,i,4)=rho(end,i)*uy(end,i)*H(end,i);
        
        %%%%% SOUTH CELLS %%%%%
        f(1,i,1)=rho(1,i)*ux(1,i);
        f(1,i,2)=rho(1,i)*ux(1,i)^2+p(1,i);
        f(1,i,3)=rho(1,i)*ux(1,i)*uy(1,i);
        f(1,i,4)=rho(1,i)*ux(1,i)*H(1,i);
        
        g(1,i,1)=rho(1,i)*uy(1,i);
        g(1,i,2)=rho(1,i)*ux(1,i)*uy(1,i);
        g(1,i,3)=rho(1,i)*uy(1,i)^2+p(1,i);
        g(1,i,4)=rho(1,i)*uy(1,i)*H(1,i);        
         
end


%%%%% CORNER CELLS %%%%%
%%%NW%%%

f(Ny-1,1,1)=rho(Ny-1,1)*ux(Ny-1,1);
f(Ny-1,1,2)=rho(Ny-1,1)*ux(Ny-1,1)^2+p(Ny-1,1);
f(Ny-1,1,3)=rho(Ny-1,1)*ux(Ny-1,1)*uy(Ny-1,1);
f(Ny-1,1,4)=rho(Ny-1,1)*ux(Ny-1,1)*H(Ny-1,1);
        
g(Ny-1,1,1)=rho(Ny-1,1)*uy(Ny-1,1);
g(Ny-1,1,2)=rho(Ny-1,1)*ux(Ny-1,1)*uy(Ny-1,1);
g(Ny-1,1,3)=rho(Ny-1,1)*uy(Ny-1,1)^2+p(Ny-1,1);
g(Ny-1,1,4)=rho(Ny-1,1)*uy(Ny-1,1)*H(Ny-1,1);

%%%NE%%%

f(Ny-1,Nx-1,1)=rho(Ny-1,Nx-1)*ux(Ny-1,Nx-1);
f(Ny-1,Nx-1,2)=rho(Ny-1,Nx-1)*ux(Ny-1,Nx-1)^2+p(Ny-1,Nx-1);
f(Ny-1,Nx-1,3)=rho(Ny-1,Nx-1)*ux(Ny-1,Nx-1)*uy(Ny-1,Nx-1);
f(Ny-1,Nx-1,4)=rho(Ny-1,Nx-1)*ux(Ny-1,Nx-1)*H(Ny-1,Nx-1);
        
g(Ny-1,Nx-1,1)=rho(Ny-1,Nx-1)*uy(Ny-1,Nx-1);
g(Ny-1,Nx-1,2)=rho(Ny-1,Nx-1)*ux(Ny-1,Nx-1)*uy(Ny-1,Nx-1);
g(Ny-1,Nx-1,3)=rho(Ny-1,Nx-1)*uy(Ny-1,Nx-1)^2+p(Ny-1,Nx-1);
g(Ny-1,Nx-1,4)=rho(Ny-1,Nx-1)*uy(Ny-1,Nx-1)*H(Ny-1,Nx-1);

%%%SE%%%%

f(1,Nx-1,1)=rho(1,Nx-1)*ux(1,Nx-1);
f(1,Nx-1,2)=rho(1,Nx-1)*ux(1,Nx-1)^2+p(1,Nx-1);
f(1,Nx-1,3)=rho(1,Nx-1)*ux(1,Nx-1)*uy(1,Nx-1);
f(1,Nx-1,4)=rho(1,Nx-1)*ux(1,Nx-1)*H(1,Nx-1);
        
g(1,Nx-1,1)=rho(1,Nx-1)*uy(1,Nx-1);
g(1,Nx-1,2)=rho(1,Nx-1)*ux(1,Nx-1)*uy(1,Nx-1);
g(1,Nx-1,3)=rho(1,Nx-1)*uy(1,Nx-1)^2+p(1,Nx-1);
g(1,Nx-1,4)=rho(1,Nx-1)*uy(1,Nx-1)*H(1,Nx-1);

%%%SW%%%%

f(1,1,1)=rho(1,1)*ux(1,1);
f(1,1,2)=rho(1,1)*ux(1,1)^2+p(1,1);
f(1,1,3)=rho(1,1)*ux(1,1)*uy(1,1);
f(1,1,4)=rho(1,1)*ux(1,1)*H(1,1);
        
g(1,1,1)=rho(1,1)*uy(1,1);
g(1,1,2)=rho(1,1)*ux(1,1)*uy(1,1);
g(1,1,3)=rho(1,1)*uy(1,1)^2+p(1,1);
g(1,1,4)=rho(1,1)*uy(1,1)*H(1,1);

end
