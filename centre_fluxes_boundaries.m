function [f,g]=centre_fluxes_boundaries(Nx,Ny,rho,u,v,H,p)

f = zeros((Ny-1),(Nx-1),4);
g = zeros((Ny-1),(Nx-1),4);


f(:,:,1)=rho.*u;
f(:,:,2)=rho.*u.^2+p;
f(:,:,3)=rho.*u.*v;
f(:,:,4)=rho.*u.*H;
        
g(:,:,1)=rho.*v;
g(:,:,2)=rho.*u.*v;
g(:,:,3)=rho.*v.^2+p;
g(:,:,4)=rho.*v.*H;

end
