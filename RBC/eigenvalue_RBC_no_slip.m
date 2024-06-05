clear all;
close all;
clc;
kx_list=2*pi/2.016; 
Ra_list=1708; %R=(n^2*pi^2+k^2)^3/k^2
Pr=1;
phi_list=0;%35/180*pi;
kappa_list=0;
% kx_list=5;
ky_list=0;

N=2;
[xi, DM] = chebdif(N+2, 4);
z=(xi+1)/2;
D2=DM(2:end-1,2:end-1,2)*2^2;
[~,dd4]=cheb4c(N+2);
D4=dd4*2^4;
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
                    inv_Laplacian=inv(Laplacian);
                    Laplacian_square=D4-2*(kx^2+ky^2)*D2+(kx^2+ky^2)^2*I;
                    
                    A=[Pr*inv_Laplacian*Laplacian_square, inv_Laplacian*Pr*Ra*(-(kx^2+ky^2));
                        I, Laplacian];
                                       
                    [eig_vec,eig_val]=eig(-A);
                    eig_val=diag(eig_val);
                    eig_val(find(real(eig_val)==Inf))=-Inf;
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

