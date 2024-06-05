clear all;
close all;
clc;
kx_list=pi/sqrt(2); 
Ra_list=27/4*pi^4+1; %R=(n^2*pi^2+k^2)^3/k^2
Pr=1;
phi_list=0;%35/180*pi;
kappa_list=0;
% kx_list=5;
ky_list=0;

N=128;
[xi, DM] = chebdif(N, 4);
z=(xi+1)/2;

D1=DM(:,:,1)*2;
D2=DM(:,:,2)*2^2;
D4=DM(:,:,4)*2^4;
I=eye(N,N);
O=zeros(N,N);
for Ra_ind=1:length(Ra_list)
    Ra=Ra_list(Ra_ind);
    for phi_ind=1:length(phi_list)
        phi=phi_list(phi_ind);
        for kappa_ind=1:length(kappa_list)
            kappa=kappa_list(kappa_ind);
            for kx_ind=1:length(kx_list)
                kx=kx_list(kx_ind);
                for ky_ind=1:length(ky_list)
                    ky=ky_list(ky_ind);
                    
                    Laplacian=D2-(kx^2+ky^2)*I;
                    Laplacian_square=D4-2*(kx^2+ky^2)*D2+(kx^2+ky^2)^2*I;
                    E=[Laplacian,O;
                        O,I];
                    A=[Pr*Laplacian_square, Pr*Ra*(-(kx^2+ky^2))*I;
                        I, Laplacian];
                    
                    %velocity boundary condition, w=0 at the boundaries
                    A(1,:)=[1,zeros(1,N-1), zeros(1,N)];
                    A(N,:)=[zeros(1,N-1), 1, zeros(1,N)];
                    E(1,:)=zeros(1,2*N);
                    E(N,:)=zeros(1,2*N);
                    
                    %velocity boundary condition, w''=0 at the boundaries
                    A(2,:)=[D2(1,:),zeros(1,N)];
                    A(N-1,:)=[D2(N,:),zeros(1,N)];
                    E(2,:)=zeros(1,2*N);
                    E(N-1,:)=zeros(1,2*N);
                    
                    %temperature boundary conditions.
                    E(1+N,:)=zeros(1,2*N);
                    E(2*N,:)=zeros(1,2*N);
                    A(1+N,:)=[zeros(1,N),[1,zeros(1,N-1)]];
                    A(2*N,:)=[zeros(1,N),zeros(1,N-1),1];

                    A_reduce=A([2:N-1,N+1:2*N-1],[2:N-1,N+1:2*N-1]);
                    E_reduce=E([2:N-1,N+1:2*N-1],[2:N-1,N+1:2*N-1]);
                    
                    [U,S,V]=svd(E);
                    
                    %[eig_vec,eig_val]=eig(A,E);
                    
                    A_svd=S^(-1/2)*U'*A*V*S^(-1/2);
                    
                    [eig_vec,eig_val]=eig(A_svd);
                    eig_val=diag(eig_val);
                    eig_val(find(real(eig_val)>10e14))=-Inf;
                    [growth_rate(Ra_ind,phi_ind,kappa_ind,kx_ind,ky_ind),max_ind]=max(real(eig_val));
                    eig_val_max{Ra_ind,phi_ind,kappa_ind,kx_ind,ky_ind}=eig_val(max_ind);
                    eig_vec_max{Ra_ind,phi_ind,kappa_ind,kx_ind,ky_ind}=eig_vec(:,max_ind);
               
                end
            end
        end
    end
end
plot(real(eig_val),imag(eig_val),'*');
xlim([-50,50]);
ylim([-50,50]);
% g_T=[1,0,0;
%     (1-kappa),kappa,0];
% g_w=[1,0,0;
%     1,0,0];
% [x_T,D2_T,D1_T,~,~]=cheb2bc(N,g_T);
% [x_w,D2_w,D1_w,~,~]=cheb2bc(N,g_w);
%  a_1 u(1)  + b_1 u'(1)  = c_1
%  a_N u(-1) + b_N u'(-1) = c_N
% g        =  boundary condition matrix = [a_1 b_1 c_1; a_N b_N c_N]

