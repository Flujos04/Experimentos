function [Nx_N,Ny_N,Nx_E,Ny_E,Nx_S,Ny_S,Nx_W,Ny_W]=face_normals(Nx,Ny,x,y)

Nx_N = zeros(Ny-1,Nx-1);
Ny_N = zeros(Ny-1,Nx-1);
Nx_E = zeros(Ny-1,Nx-1);
Ny_E = zeros(Ny-1,Nx-1);
Nx_S = zeros(Ny-1,Nx-1);
Ny_S = zeros(Ny-1,Nx-1);
Nx_W = zeros(Ny-1,Nx-1);
Ny_W = zeros(Ny-1,Nx-1);

for i=1:(Nx-1)
    for j=1:(Ny-1)
        
        %North
        Nx_N(j,i)=y(j+1,i)-y(j+1,i+1);
        Ny_N(j,i)=-(x(i)-x(i+1));
        
        %East
        Nx_E(j,i)=y(j+1,i+1)-y(j,i+1);
        Ny_E(j,i)=-(x(i+1)-x(i+1));
        
        %South
        Nx_S(j,i)=y(j,i+1)-y(j,i);
        Ny_S(j,i)=-(x(i+1)-x(i));
        
        %West
        Nx_W(j,i)=y(j,i)-y(j+1,i);
        Ny_W(j,i)=-(x(i)-x(i));
    end
end
