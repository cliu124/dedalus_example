clear all;
close all;
clc;

%Only change this parameter
R_rho_list=0.1;

tau_list=0.01;
Wst_list=0;
Pr_list=7;

k_list=linspace(-3,3,201);
l_list=logspace(-7,1,401);
phi_list=1;

%This is reproduction for validation
%grad_T_vertical=1;
%phi_C=phi_list;

%This corresponds to horizontal gradient only
grad_T_vertical=0;
phi_C=0;

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
                            phi*k/l-grad_T_vertical*1, -K2, 0;
                            phi_C*k/l-1/R_rho,0,-tau*K2+Wst*1i*k];
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