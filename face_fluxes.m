function [FN, FE, FS, FW] = face_fluxes(f,g,Nx_N,Ny_N,Nx_E,Ny_E,Nx_S,Ny_S,Nx_W,Ny_W,Ny,Nx,D_north,D_west,D_east,D_south)


%NORTH FACES
FN = zeros((Ny-1),(Nx-1),4);
    %Boundary
        j=Ny-1;
        FN(j,:,:)=f(j,:,:).*Nx_N(j,:)+g(j,:,:).*Ny_N(j,:)-D_north(j,:,:);
                
    %Rest        
        for i=1:(Nx-1)
          for j=1:(Ny-2)
               FN (j,i,:) =(1/2)*(f(j,i,:)+f(j+1,i,:))*Nx_N(j,i) + (1/2)*(g(j,i,:)+g(j+1,i,:))*Ny_N(j,i)-D_north(j,i,:); 
          end
        end
        
        
%EAST FACES
FE = zeros((Ny-1),(Nx-1),4);
    %Boundary
        i=Nx-1;
        FE(:,i,:)=f(:,i,:).*Nx_E(:,i)+g(:,i,:).*Ny_E(:,i)-D_east(:,i,:);
                
    %Rest        
        for i=1:(Nx-2)
          for j=1:(Ny-1)
               FE (j,i,:) =(1/2)*(f(j,i,:)+f(j,i+1,:))*Nx_E(j,i) + (1/2)*(g(j,i,:)+g(j,i+1,:))*Ny_E(j,i)-D_east(j,i,:); 
          end
        end
        
        
%SOUTH FACES
FS = zeros((Ny-1),(Nx-1),4);
    %Boundary
        j=1;
        FN(j,:,:)=f(j,:,:).*Nx_S(j,:)+g(j,:,:).*Ny_S(j,:)-D_south(j,:,:);
                
    %Rest        
        for i=1:(Nx-1)
          for j=2:(Ny-1)
               FS (j,i,:) =(1/2)*(f(j,i,:)+f(j-1,i,:))*Nx_S(j,i) + (1/2)*(g(j,i,:)+g(j-1,i,:))*Ny_S(j,i)-D_south(j,i,:);
          end
        end        
 
        
%WEST FACES
FW = zeros((Ny-1),(Nx-1),4);
    %Boundary
        i=1;
        FW(:,i,:)=f(:,i,:).*Nx_W(:,i)+g(:,i,:).*Ny_W(:,i)-D_west(:,i,:);
                
    %Rest        
        for i=2:(Nx-1)
          for j=1:(Ny-1)
               FW (j,i,:) =(1/2)*(f(j,i,:)+f(j,i-1,:))*Nx_W(j,i) + (1/2)*(g(j,i,:)+g(j,i-1,:))*Ny_W(j,i)-D_west(j,i,:);
          end
        end 
end
