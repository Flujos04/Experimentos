function [f,g]=centre_fluxes_boundaries(Nx,Ny,rho,ux,uy,H,p,rho_w,u_w,p_w,v_w,H_w,rho_e,u_e,p_e,v_e,H_e,v_n,v_s)


%We are going to work with f and g in matrixes. Then, they will be converted to vectors for simplicity
f = zeros((Ny-1),(Nx-1),4);
g = zeros((Ny-1),(Nx-1),4);


% f(:,:,1)=rho.*ux;
% f(:,:,2)=rho.*ux.^2+p;
% f(:,:,3)=rho.*ux.*uy;
% f(:,:,4)=rho.*ux.*H;
%         
% g(:,:,1)=rho.*uy;
% g(:,:,2)=rho.*ux.*uy;
% g(:,:,3)=rho.*uy.^2+p;
% g(:,:,4)=rho.*uy.*H;



% f(:,:,1)=U(:,:,2);
% f(:,:,2)=(U(:,:,2).^2)./U(:,:,1)+p(:,:);
% f(:,:,3)=U(:,:,2).*U(:,:,3)./U(:,:,1);
% f(:,:,4)=U(:,:,2).*H(:,:);
%         
% g(:,:,1)=U(:,:,3);
% g(:,:,2)=U(:,:,2).*U(:,:,3)./U(:,:,1);
% g(:,:,3)=(U(:,:,3).^2)./U(:,:,1)+p(:,:);
% g(:,:,4)=U(:,:,3).*H(:,:);



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
        f(j,1,1)=rho_w(j,1)*u_w(j,1);
        f(j,1,2)=rho_w(j,1)*u_w(j,1)^2+p_w(j,1);
        f(j,1,3)=rho_w(j,1)*u_w(j,1)*v_w(j,1);
        f(j,1,4)=rho_w(j,1)*u_w(j,1)*H_w(j,1);
        
        g(j,1,1)=rho_w(j,1)*v_w(j,1);
        g(j,1,2)=rho_w(j,1)*u_w(j,1)*v_w(j,1);
        g(j,1,3)=rho_w(j,1)*v_w(j,1)^2+p_w(j,1);
        g(j,1,4)=rho_w(j,1)*v_w(j,1)*H_w(j,1);
        
        %%%%% EAST CELLS %%%%%
        f(j,end,1)=rho_e(j,1)*u_e(j,1);
        f(j,end,2)=rho_e(j,1)*u_e(j,1)^2+p_e(j,1);
        f(j,end,3)=rho_e(j,1)*u_e(j,1)*v_e(j,1);
        f(j,end,4)=rho_e(j,1)*u_e(j,1)*H_e(j,1);
        
        g(j,end,1)=rho_e(j,1)*v_e(j,1);
        g(j,end,2)=rho_e(j,1)*u_e(j,1)*v_e(j,1);
        g(j,end,3)=rho_e(j,1)*v_e(j,1)^2+p_e(j,1);
        g(j,end,4)=rho_e(j,1)*v_e(j,1)*H_e(j,1);
        
    end
    
        %%%%% NORTH CELLS %%%%%
        f(end,i,1)=rho(end,i)*ux(end,i);
        f(end,i,2)=rho(end,i)*ux(end,i)^2+p(end,i);
        f(end,i,3)=rho(end,i)*ux(end,i)*v_n(i,1);
        f(end,i,4)=rho(end,i)*ux(end,i)*H(end,i);
        
        g(end,i,1)=rho(end,i)*v_n(i,1);
        g(end,i,2)=rho(end,i)*ux(end,i)*v_n(i,1);
        g(end,i,3)=rho(end,i)*v_n(i,1)^2+p(end,i);
        g(end,i,4)=rho(end,i)*v_n(i,1)*H(end,i);
        
        %%%%% SOUTH CELLS %%%%%
        f(1,i,1)=rho(1,i)*ux(1,i);
        f(1,i,2)=rho(1,i)*ux(1,i)^2+p(1,i);
        f(1,i,3)=rho(1,i)*ux(1,i)*v_s(i,1);
        f(1,i,4)=rho(1,i)*ux(1,i)*H(1,i);
        
        g(1,i,1)=rho(1,i)*v_s(i,1);
        g(1,i,2)=rho(1,i)*ux(1,i)*v_s(i,1);
        g(1,i,3)=rho(1,i)*v_s(i,1)^2+p(1,i);
        g(1,i,4)=rho(1,i)*v_s(i,1)*H(1,i);        
         
end


%%%%% CORNER CELLS %%%%%
%%%NW%%%

f(Ny-1,1,1)=rho_w(Ny-1,1)*u_w(Ny-1,1);
f(Ny-1,1,2)=rho_w(Ny-1,1)*u_w(Ny-1,1)^2+p_w(Ny-1,1);
f(Ny-1,1,3)=rho_w(Ny-1,1)*u_w(Ny-1,1)*v_w(Ny-1,1);
f(Ny-1,1,4)=rho_w(Ny-1,1)*u_w(Ny-1,1)*H_w(Ny-1,1);
        
g(Ny-1,1,1)=rho_w(Ny-1,1)*v_w(Ny-1,1);
g(Ny-1,1,2)=rho_w(Ny-1,1)*u_w(Ny-1,1)*v_w(Ny-1,1);
g(Ny-1,1,3)=rho_w(Ny-1,1)*v_w(Ny-1,1)^2+p_w(Ny-1,1);
g(Ny-1,1,4)=rho_w(Ny-1,1)*v_w(Ny-1,1)*H_w(Ny-1,1);

%%%NE%%%

f(Ny-1,Nx-1,1)=rho_e(Ny-1,1)*u_e(Ny-1,1);
f(Ny-1,Nx-1,2)=rho_e(Ny-1,1)*u_e(Ny-1,1)^2+p_e(Ny-1,1);
f(Ny-1,Nx-1,3)=rho_e(Ny-1,1)*u_e(Ny-1,1)*v_e(Ny-1,1);
f(Ny-1,Nx-1,4)=rho_e(Ny-1,1)*u_e(Ny-1,1)*H_e(Ny-1,1);
        
g(Ny-1,Nx-1,1)=rho_e(Ny-1,1)*v_e(Ny-1,1);
g(Ny-1,Nx-1,2)=rho_e(Ny-1,1)*u_e(Ny-1,1)*v_e(Ny-1,1);
g(Ny-1,Nx-1,3)=rho_e(Ny-1,1)*v_e(Ny-1,1)^2+p_e(Ny-1,1);
g(Ny-1,Nx-1,4)=rho_e(Ny-1,1)*v_e(Ny-1,1)*H_e(Ny-1,1);


%%%SE%%%%

f(1,Nx-1,1)=rho_e(1,1)*u_e(1,1);
f(1,Nx-1,2)=rho_e(1,1)*u_e(1,1)^2+p_e(1,1);
f(1,Nx-1,3)=rho_e(1,1)*u_e(1,1)*v_s(1,1);
f(1,Nx-1,4)=rho_e(1,1)*u_e(1,1)*H_e(1,1);
        
g(1,Nx-1,1)=rho_e(1,1)*v_s(1,1);
g(1,Nx-1,2)=rho_e(1,1)*u_e(1,1)*v_s(1,1);
g(1,Nx-1,3)=rho_e(1,1)*v_s(1,1)^2+p_e(1,1);
g(1,Nx-1,4)=rho_e(1,1)*v_s(1,1)*H_e(1,1);


%%%SW%%%%

f(1,1,1)=rho_w(1,1)*u_w(1,1);
f(1,1,2)=rho_w(1,1)*u_w(1,1)^2+p_w(1,1);
f(1,1,3)=rho_w(1,1)*u_w(1,1)*v_w(1,1);
f(1,1,4)=rho_w(1,1)*u_w(1,1)*H_w(1,1);
        
g(1,1,1)=rho_w(1,1)*v_w(1,1);
g(1,1,2)=rho_w(1,1)*u_w(1,1)*v_w(1,1);
g(1,1,3)=rho_w(1,1)*v_w(1,1)^2+p_w(1,1);
g(1,1,4)=rho_w(1,1)*v_w(1,1)*H_w(1,1);



end
