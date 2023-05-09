clear all;
close all;
clc;

k_list=linspace(-1,1,101);
l_list=logspace(-5,1,201);
phi_list=0;
R_rho_list=0.001;%;[0.01,0.1,0.5,1,10,50,100];
tau_list=278.767;
Wst_list=1456;
Pr_list=700;

for Pr_ind=1:length(Pr_list)
    Pr=Pr_list(Pr_ind);
   for Wst_ind=1:length(Wst_list)
       Wst=Wst_list(Wst_ind);
      for tau_ind=1:length(tau_list)
          tau=tau_list(tau_ind);
          for R_rho_ind=1:length(R_rho_list)
             R_rho=R_rho_list(R_rho_ind);
             for phi_ind=1:length(phi_list)
                 phi=phi_list(phi_ind);
                 for k_ind=1:length(k_list)
                    k=k_list(k_ind);
                    for l_ind=1:length(l_list)
                        l=l_list(l_ind);
                        K=sqrt(k^2+l^2);
                        K2=k^2+l^2;
                        M=[K2/l^2/Pr, 0,0;
                            0,1,0;
                            0,0,1];
                        A=[-K2*K2/l^2, 1,-1;
                            phi*k/l-1, -K2, 0;
                            phi*k/l-1/R_rho,0,-tau*K2+Wst*1i*k];
                        [eig_vec,eig_val]=eig(A,M);
                        growth_rate{Pr_ind,Wst_ind,tau_ind,R_rho_ind,phi_ind}(k_ind,l_ind)=...
                            max(real(diag(eig_val)));
                    end
                 end
             end
          end
      end
   end
end

growth_rate{1}(find(growth_rate{1}<0))=NaN;
data{1}.x=k_list;
data{1}.y=log10(l_list);
data{1}.z=growth_rate{1}';
plot_config.label_list={1,'$k$','log$_{10}(l)$'};
plot_contour(data,plot_config);