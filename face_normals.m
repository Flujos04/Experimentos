function [dS]=face_normals(Nx,Ny,x,y)

dS=cell(Nx-1,Ny-1);

for i=1:(Nx-1)
    for j=1:(Ny-1)
        dS{i,j}=zeros(4,2);
    end
end

for i=1:(Nx-1)
    for j=1:(Ny-1)
        
        %North
        dS{i,j}(1,1)=y(j+1,i)-y(j+1,i+1);
        dS{i,j}(1,2)=-(x(i)-x(i+1));
        
        %East
        dS{i,j}(2,1)=y(j+1,i+1)-y(j,i+1);
        dS{i,j}(2,2)=-(x(i+1)-x(i+1));
        
        %South
        dS{i,j}(3,1)=y(j,i+1)-y(i,j);
        dS{i,j}(3,2)=-(x(i+1)-x(i));
        
        %West
        dS{i,j}(4,1)=y(i,j)-y(j+1,i);
        dS{i,j}(4,2)=-(x(i)-x(i));
    end
end
