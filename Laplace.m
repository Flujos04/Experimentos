clc;clear;close all;

dx=0.01;
dy=0.01;
L=1; 

xaxis=0:dx:L; 
yaxis=0:dy:L; 

Nx=L/dx + 1;
Ny=L/dy + 1;

xx=zeros(Nx^2,1);
A=zeros(Nx^2,Ny^2);
b=zeros(Nx^2,1);
U=zeros(Nx^2);



%%%%%%%%%% 5STENCIL %%%%%%%%%%%%%%

for i=2:(Nx-1)
    for j=2:(Ny-1)
        k=i+Nx*(j-1);
        Arow=zeros(1,Nx^2);
        Arow(k)=-4;
        Arow(k-1)=1;
        Arow(k+1)=1;
        Arow(k-Nx)=1;
        Arow(k+Nx)=1;
        A(k,:)=Arow;
    end
end

%%%%%%%% SOUTH-BOUNDARY %%%%%%%%%%%%%%%

for j=1
    for i=1:Nx
        k=i+Nx*(j-1);
        Arow=zeros(1,Nx^2);
        Arow(k)=1;
        b(k)=0;
        A(k,:)=Arow;
    end
end

%%%%%%%%% EASTBOUNDARY %%%%%%%%%%%%%%

for i=Nx
    for j=2:(Nx-1)
        k=i+Nx*(j-1);
        Arow=zeros(1,Nx^2);
        Arow(k)=1;
        b(k)=0;
        A(k,:)=Arow;
    end
end

%%%%%%%%% WEST BOUNDARY %%%%%%%%%%%%


for i=1
    for j=2:(Nx-1)
        k=i+Nx*(j-1);
        Arow=zeros(1,Nx^2);
        Arow(k)=1;
        b(k)=0;
        A(k,:)=Arow;
    end
end


%%%%%%%% NORTH BOUNDARY %%%%%%%%%%%

funcnorth1=@(xx) sin(pi*xx);
fnorth=funcnorth1(xaxis);

for j=Ny
    for i=1:Nx
        k=i+Nx*(j-1);
        Arow=zeros(1,Nx^2);
        Arow(k)=1;
        b(k)=fnorth(i);
        A(k,:)=Arow;
        pene(i)=sin(pi*xaxis(i));
    end
end


xx=mldivide(A,b);
solution=zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
         k=i+Nx*(j-1);
         solution(i,j)=xx(k);
    end
end

for i=1:Nx
    for j=1:Ny
        analytical(i,j)=(sinh(pi.*yaxis(j)).*sin(pi.*xaxis(i)))./(sinh(pi));
    end
end

figure; 
surf(yaxis,xaxis,solution); 
title('Discrete solution'); shading flat; axis equal;
xlabel('y');
ylabel('x');
zlabel('function');

figure; 
surf(yaxis,xaxis,analytical); 
title('Analytical Solution'); shading flat; axis equal;
xlabel('y');
ylabel('x');
zlabel('function');
     





