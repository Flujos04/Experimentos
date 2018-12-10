
function [x,y,y0,dx,dy]=mesher(L,Nx,Ny)


dx=3*L/(Nx-1);
x=0:dx:3*L;
R=(1.3)*L;
xarc=L:dx:L*2;
yarc=sqrt(R^2 - (xarc-1.5*L).^2) - 1.2*L;
    
y0=[(0:dx:L-dx)*0 yarc (0:dx:L-dx)*0];
dy=(L-y0)./(Ny-1);
y=zeros(Ny,Nx);
    
for i=1:Nx
    y(:,i)=[y0(i):dy(i):L]';
end

for i=1:Ny
    
   xmatrix(i,:)=x; 
end
figure;hold on;
surf(xmatrix,y,xmatrix.*0);
end
