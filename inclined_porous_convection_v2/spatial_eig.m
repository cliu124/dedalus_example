clear all;
close all;
clc;

Ra_list=100;%linspace(25,30,20);
phi_list=35/180*pi;
kappa_list=0;
lambda_t=0;

N=128;
[xi, DM] = chebdif(N, 2);
z=(xi+1)/2;
D1=DM(:,:,1)*2;
D2=DM(:,:,2)*2^2;
I=eye(N,N);
O=zeros(N,N);
for Ra_ind=1:length(Ra_list)
    Ra=Ra_list(Ra_ind);
    for phi_ind=1:length(phi_list)
        phi=phi_list(phi_ind);
        for kappa_ind=1:length(kappa_list)
            kappa=kappa_list(kappa_ind);
                    
            E=blkdiag(I,I,I,I);
            A=[O, I, O, O;
                -D2, O, Ra*sin(phi)*D1, -Ra*cos(phi)*I;
                O, O, O, I;
                O, I, lambda_t-D2, Ra*sin(phi)*diag(1/2-z)];

            %velocity boundary condition, w(z=0)=w(z=1)=0
            A(1+N,:)=[zeros(1,N), 1,zeros(1,N-1), zeros(1,N), zeros(1,N)];
            A(2*N,:)=[zeros(1,N), zeros(1,N-1),1, zeros(1,N), zeros(1,N)];
            E(1+N,:)=zeros(1,4*N);
            E(2*N,:)=zeros(1,4*N);
            %temperature boundary conditions.
            E(1+2*N,:)=zeros(1,4*N);
            E(3*N,:)=zeros(1,4*N);
            A(1+2*N,:)=[zeros(1,N),zeros(1,N),((1-kappa)*[1,zeros(1,N-1)]+kappa*D1(1,:)),zeros(1,N)];
            A(3*N,:)=[zeros(1,N),zeros(1,N),zeros(1,N-1),1,zeros(1,N)];

            [eig_vec,eig_val]=eig(A,E);
            eig_val=diag(eig_val);
            [~,sort_index]=sort(abs(real(eig_val)));
            eig_sort=eig_val(sort_index);
            eig_sort(find(eig_sort==0))=[];
        end
    end
end
N_disp=10;
plot(real(eig_sort),imag(eig_sort),'*');
xlim([min(real(eig_sort(1:N_disp))),max(real(eig_sort(1:N_disp)))]);
ylim([min(imag(eig_sort(1:N_disp))),max(imag(eig_sort(1:N_disp)))]);
disp(['First ',num2str(N_disp),' spatial eigenvalue with real part cloest to the origin:']);
eig_sort(1:N_disp);
