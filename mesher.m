
function [x,y,y0,dx,dy]=mesher(L,Nx,Ny)


dx=3*L/(Nx-1);
x=0:dx:3*L;
R=(1.3)*L;

xinit=0:dx:L;
xarc=L:dx:L*2;
xfinal=2*L:dx:3*L;

x(1:length(xinit))=xinit; 
x(length(xinit):length(xinot)+length(xarc)-1)=xarc;
x(length(xint)+length(xarc)-1:length(xinit)+length(xarc)+length(xfinal)-2)=xfinal;



    lilbump=sqrt(R^2 - (xarc-1.5*L).^2) - 1.2*L;
    y0(length(xinit):length(xinit)+length(xarc)-1)=lilbump;
    y0(1:length(xinit))=0;
    y0(length(xinit)+length(xarc)-1:length(xinit)+length(xarc)+length(xfinal)-2)=0;
    dy=(L-y0)./(Ny-1);
    y=zeros(Ny,Nx);
    
for i=1:Nx
    y(:,i)=[y0(i):dy(i):L]';
end

for i=1:Ny
    
   xmatrix(i,:)=x; 
end

surf(xmatrix,y,xmatrix.*0);
end
