clear all;
close all;
clc;

Ra_list=100;%linspace(25,30,20);
phi_list=35/180*pi;
kappa_list=0:0.001:0.01;
c=0;

N=32;
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
                   
            %Update c here if you need to change c as a function of kappa
            %c=kappa;
            clear E A; 
            E=blkdiag(I,I,I,I);
            A=[O, I, O, O;
                -D2, O, Ra*sin(phi)*D1, -Ra*cos(phi)*I;
                O, O, O, I;
                O, I, -D2, Ra*sin(phi)*diag(1/2-z)-c*I];

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
            
            %Update the eigenshuffle to track the change over kappa
            %Update 2023/11/05
            Aseq(:,:,kappa_ind)=A;
            Eseq(:,:,kappa_ind)=E;
        end
    end
end
N_disp=20;
plot(real(eig_sort),imag(eig_sort),'*');
xlim([min(real(eig_sort(1:N_disp))),max(real(eig_sort(1:N_disp)))]);
ylim([min(imag(eig_sort(1:N_disp))),max(imag(eig_sort(1:N_disp)))]);
disp(['First ',num2str(N_disp),' spatial eigenvalue with real part cloest to the origin:']);
eig_sort(1:N_disp)


%Add eigenshuffle to track the change of eigenvalues, 2023/11/05
[Vseq,Dseq]=eigenshuffle(Aseq,Eseq);
[~,sort_index]=sort(abs(real(Dseq(:,1))));
Dseq_sort=Dseq(sort_index,:);
ind=find(Dseq_sort(:,1)==0);
Dseq_sort(ind,:)=[];

figure(1)
plot(kappa_list,real(Dseq_sort(1,:))); hold on;
plot(kappa_list,imag(Dseq_sort(1,:)));

figure(2)
plot(kappa_list,real(Dseq_sort(2,:))); hold on;
plot(kappa_list,imag(Dseq_sort(2,:)));

figure(3)
plot(kappa_list,real(Dseq_sort(3,:))); hold on;
plot(kappa_list,imag(Dseq_sort(3,:)));

figure(4)
plot(kappa_list,real(Dseq_sort(4,:))); hold on;
plot(kappa_list,imag(Dseq_sort(4,:)));
