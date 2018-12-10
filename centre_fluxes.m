function [f,g]=centre_fluxes(Nx,Ny,rho,ux,uy,H,p)

%We are going to work with F and G in matrixes. Then, they will be converted to vectors for simplicity
f = zeros((Ny-1),(Nx-1),4);
g = zeros((Ny-1),(Nx-1),4);

for i=1:(Nx-1)
    for j=1:(Ny-1)
        f(j,i,1)=rho(j,i)*ux(j,i);
        f(j,i,2)=rho(j,i)*ux(j,i)^2+p(j,i);
        f(j,i,3)=rho(j,i)*ux(j,i)*uy(j,i);
        f(j,i,4)=rho(j,i)*ux(j,i)*H(j,i);
        
        g(j,i,1)=rho(j,i)*uy(j,i);
        g(j,i,2)=rho(j,i)*ux(j,i)*uy(j,i);
        g(j,i,3)=rho(j,i)*uy(j,i)^2+p(j,i);
        g(j,i,4)=rho(j,i)*uy(j,i)*H(j,i);
    end
end

end
