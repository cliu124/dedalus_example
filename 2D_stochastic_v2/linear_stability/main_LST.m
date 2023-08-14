clear all;
close all;
clc;

flag.name='growth_rate';
flag.no_flux_com=0;
flag.no_shear_com=0;
flag.solve='eig'; %'finished' if just want to load data and plot
flag.viscous_unit=0;
switch flag.name
    case 'growth_rate'
        n_elevator_list=1;
%         Ra_T_q_list=8*10^8;
        %Ra_T_q_list=46892.1;
        %Pr_list=1;

        Ra_T_q_list=46892;
        Pr_list=1;
        Lx_list=0.1*2*pi;
%         kz_list=(0.01:0.5:2*n_elevator)*2*pi;
%         kz_list=(0.01:0.1:2)*2*pi;
        kz_list=[1]*2*pi;
%         kz_list=2*pi;
        flag.no_flux_com=0; %compare the growth rate without flux feedback
        flag.no_shear_com=0; %compare the growth rate without shear flow
        flag.viscous_unit=0;
        flag.plot_spectrum=1;
    case 'Lx'
        n_elevator_list=1;
        Ra_T_q_list=10^8;
        Pr_list=1;
        Lx_list=[0.05,0.1:0.1:1]*2*pi;
        kz_list=(0.01:0.01:2)*2*pi;
    case 'Ra_T_q'
        n_elevator_list=1;
        Ra_T_q_list=[20000,30000,40000,60000,10^5,10^6,10^7,10^8];
        Pr_list=1;
        Lx_list=2*pi/10;
        kz_list=(0.01:0.01:2)*2*pi;
    case 'Pr'
        n_elevator_list=1;
        Ra_T_q_list=10^8;
        Pr_list=[0.01,0.1,1,7,70,700];
        Lx_list=2*pi/10;
        kz_list=(0.01:0.01:2)*2*pi;
    case 'ke'
        n_elevator_list=[1,2,3,4,5,6];
        Ra_T_q_list=10^8;
        Pr_list=1;
        Lx_list=2*pi/5;
        kz_list=(0.01:0.01:5)*2*pi;
end

%kz_list_int=0.01;
%kz_list_num=0.01:0.02:10;%logspace(-3,1,30);

% kz_list_int=0.001*2*pi;%logspace(-6,-3,4)*(2*pi);
% kz_list_int=10^(-12);
% kz_list_num=logspace(-4,0,5);
% kz_list=[kz_list_int,kz_list_num];

Nx=128;
ky_list=0;
%get the differential matrix. 
[x_2pi,D1x_2pi] = fourdif(Nx,1);
[~,D2x_2pi] = fourdif(Nx,2);
[~,D4x_2pi] = fourdif(Nx,4);

% D1x_2pi_fft=fourdifft(x_2pi,1);
% D2x_2pi_fft=fourdifft(x_2pi,2);

I=eye(Nx,Nx);
O=zeros(Nx,Nx);
flux_T=1;
if ~strcmp(flag.solve,'finished')
    for Lx_ind=1:length(Lx_list)
        Lx=Lx_list(Lx_ind);
        for n_elevator_ind=1:length(n_elevator_list)
            n_elevator=n_elevator_list(n_elevator_ind);
            kx_mean=n_elevator*2*pi/Lx;
            x=x_2pi/2/pi*Lx;
            D1x=D1x_2pi*2*pi/Lx;
            D2x=D2x_2pi*(2*pi/Lx)^2;
            D4x=D4x_2pi*(2*pi/Lx)^4;
            dx=diff(x); dx=dx(1);
            for Ra_T_q_ind=1:length(Ra_T_q_list)
                Ra_T_q=Ra_T_q_list(Ra_T_q_ind);
                w_hat=sqrt(-kx_mean^2/2+Ra_T_q/2/kx_mean^2);
                T_hat=kx_mean^2/Ra_T_q*w_hat;
                %dy_T_mean_q=kx_mean^4/Ra_T_q;
                W_mean=2*w_hat*cos(kx_mean*x);
                d_W_mean=-2*w_hat*kx_mean*sin(kx_mean*x);
                dd_W_mean=-2*w_hat*kx_mean^2*cos(kx_mean*x);
                %T_mean=2*T_hat*cos(kx_mean*x);
                d_T_mean=-2*T_hat*kx_mean*sin(kx_mean*x);
                wT_int=1-kx_mean^4/Ra_T_q;
                for Pr_ind=1:length(Pr_list)
                    Pr=Pr_list(Pr_ind);
                    for ky_ind=1:length(ky_list)
                        ky=ky_list(ky_ind);
                        for kz_ind=1:length(kz_list)
                            kz=kz_list(kz_ind);
                            Laplacian=D2x-(ky^2+kz^2)*I;
                            flag_kz_0=0;%(kz==0);

                            %below is the formulation using velocity vorticity ones
                            Laplacian_square=D4x-2*(ky^2+kz^2)*D2x+(ky^2+kz^2)^2*I;
                            Laplacian_inv=pinv(Laplacian);
                            A11=Laplacian_inv*(-1i*kz*diag(W_mean)*Laplacian+1i*kz*diag(dd_W_mean)+Pr*Laplacian_square);
                            A12=O;
                            A13=Pr*Ra_T_q*Laplacian_inv*(-1i*kz*D1x);

                            A21=-1i*ky*diag(d_W_mean);
                            A22=-1i*kz*diag(W_mean)+Pr*Laplacian;
                            A23=Pr*Ra_T_q*1i*ky*I;

                            A31=-diag(d_T_mean)+(1-flux_T*wT_int)*1i*kz*D1x/(ky^2+kz^2);
                            A32=(1-flux_T*wT_int)*(-1i*ky)/(ky^2+kz^2)*I;
                            A33=-1i*kz*diag(W_mean)+Laplacian;

                            A=[A11,A12,A13;
                                A21,A22,A23;
                                A31,A32,A33];

                            M=blkdiag(I,I,I);

                            if flag.viscous_unit
                                if Pr~=0
                                    w_hat=sqrt(-kx_mean^2/2+Ra_T_q/2/kx_mean^2)/Pr;
                                    wT_int=(1-kx_mean^4/Ra_T_q)/Pr^2;
                                else
                                    %Pr=0, this is singular limit.
                                    w_hat=sqrt(-kx_mean^2/2+Ra_T_q/2/kx_mean^2);
                                    wT_int=(1-kx_mean^4/Ra_T_q);
                                end
                                T_hat=kx_mean^2/Ra_T_q*w_hat; %Note that this is computed based on w_hat, so no need to divied by Pr again
                                %dy_T_mean_q=kx_mean^4/Ra_T_q;
                                W_mean=2*w_hat*cos(kx_mean*x);
                                d_W_mean=-2*w_hat*kx_mean*sin(kx_mean*x);
                                dd_W_mean=-2*w_hat*kx_mean^2*cos(kx_mean*x);
                                %T_mean=2*T_hat*cos(kx_mean*x);
                                d_T_mean=-2*T_hat*kx_mean*sin(kx_mean*x);

                                A11=Laplacian_inv*(-1i*kz*diag(W_mean)*Laplacian+1i*kz*diag(dd_W_mean)+Laplacian_square);
                                A12=O;
                                A13=Ra_T_q*Laplacian_inv*(-1i*kz*D1x);

                                A21=-1i*ky*diag(d_W_mean);
                                A22=-1i*kz*diag(W_mean)+Laplacian;
                                A23=Ra_T_q*1i*ky*I;

                                A31=-diag(d_T_mean*Pr)+(1-flux_T*wT_int*Pr^2)*1i*kz*D1x/(ky^2+kz^2);
                                A32=(1-flux_T*wT_int*Pr^2)*(-1i*ky)/(ky^2+kz^2)*I;
                                A33=-1i*kz*diag(W_mean*Pr)+Laplacian;

                                A=[A11,A12,A13;
                                    A21,A22,A23;
                                    A31,A32,A33];

                                M=blkdiag(I,I,Pr*I);
                            end


                        %corresponding to u, v, w, p T, dy_T_mean_q
                    %     A=[-1i*kz*diag(W_mean)+Pr*Laplacian,O,O,-D1x,O;
                    %             O,-1i*kz*diag(W_mean)+Pr*Laplacian,O,-1i*ky*I,O;
                    %             -diag(d_W_mean),O,-1i*kz*diag(W_mean)+Pr*Laplacian,-1i*kz*I,Pr*Ra_T_q*I;
                    %             D1x, 1i*ky*I, 1i*kz*I, O, O;
                    %             -diag(d_T_mean),O, (1-flux_T*wT_int)*I-flux_T*flag_kz_0*diag(W_mean)*ones(Nx,1)*T_mean', O, -1i*kz*diag(W_mean)+Laplacian-flux_T*flag_kz_0*diag(W_mean)*ones(Nx,1)*W_mean'];
                    %     M=blkdiag(I,I,I,O,I);


                    %     if kz<0.01
                    %        A=[-1i*kz*diag(W_mean)+Pr*Laplacian,O,O,zeros(Nx,1);
                    %             O,-1i*kz*diag(W_mean)+Pr*Laplacian,Pr*Ra_T_q*I,zeros(Nx,1);
                    %             O, (1-flux_T*wT_int)*I-flux_T*flag_kz_0*diag(W_mean)*ones(Nx,1)*T_mean', -1i*kz*diag(W_mean)+Laplacian-flux_T*flag_kz_0*diag(W_mean)*ones(Nx,1)*W_mean',ones(Nx,1);
                    %             zeros(1,Nx),zeros(1,Nx),ones(1,Nx),0];
                    %        M=blkdiag(I,I,I,0);
                    %     end

                    %     A=[-1i*kz*diag(W_mean)+Pr*Laplacian,O,O,-D1x,O,zeros(Nx,1);
                    %             O,-1i*kz*diag(W_mean)+Pr*Laplacian,O,-1i*ky*I,O,zeros(Nx,1);
                    %             -diag(d_W_mean),O,-1i*kz*diag(W_mean)+Pr*Laplacian,-1i*kz*I,Pr*Ra_T_q*I,zeros(Nx,1);
                    %             D1x, 1i*ky*I, 1i*kz*I, O, O,sawtooth'*dx;
                    %             -diag(d_T_mean),O, (1-flux_T*wT_int)*I-flux_T*flag_kz_0*diag(W_mean)*ones(Nx,1)*T_mean', O, -1i*kz*diag(W_mean)+Laplacian-flux_T*flag_kz_0*diag(W_mean)*ones(Nx,1)*W_mean',zeros(Nx,1);
                    %             zeros(1,Nx),zeros(1,Nx),zeros(1,Nx),sawtooth*dx,zeros(1,Nx),0];

                    %         A=[];% u
                                 % omega_x
                                 %T
                            [eig_vec,eig_val]=eig(A,M);
                            eig_val=diag(eig_val);
                            eig_val(find(real(eig_val)==Inf))=-Inf;
                            [growth_rate(Lx_ind,n_elevator_ind,Ra_T_q_ind,Pr_ind,ky_ind,kz_ind),max_ind]=max(real(eig_val));
                            eig_val_max{Lx_ind,n_elevator_ind,Ra_T_q_ind,Pr_ind,ky_ind,kz_ind}=eig_val(max_ind);
                            eig_vec_max{Lx_ind,n_elevator_ind,Ra_T_q_ind,Pr_ind,ky_ind,kz_ind}=eig_vec(:,max_ind);

                            if flag.plot_spectrum
                                data{1}.x=real(eig_val);
                                data{1}.y=imag(eig_val);
                                data{2}.x=0;
                                data{2}.y=0;
                                plot_config.xlim_list=[1,-200,50];
                                plot_config.ylim_list=[1,-150,150];
                                plot_config.Markerindex=3;
                                plot_config.user_color_style_marker_list={'k*','bo'};
                                plot_config.label_list={1,'Re($\lambda)$','Im($\lambda$)'};
                                plot_config.name=['spectrum',flag.name,'.png'];
                                plot_line(data,plot_config);
                            end
                            
                            if flag.no_flux_com
                                flux_T=0;

                                A31=-diag(d_T_mean)+(1-flux_T*wT_int)*1i*kz*D1x/(ky^2+kz^2);
                                A32=(1-flux_T*wT_int)*(-1i*ky)/(ky^2+kz^2)*I;

                                A=[A11,A12,A13;
                                    A21,A22,A23;
                                    A31,A32,A33];

                                [eig_vec,eig_val]=eig(A,M);
                                eig_val=diag(eig_val);
                                eig_val(find(real(eig_val)==Inf))=-Inf;
                                [growth_rate_no_flux_com(Lx_ind,n_elevator_ind,Ra_T_q_ind,Pr_ind,ky_ind,kz_ind),max_ind]=max(real(eig_val));

                            end

                            if flag.no_shear_com
                                shear=0;
                                flux_T=1;
                                A11=Laplacian_inv*(-1i*kz*diag(W_mean)*shear*Laplacian+1i*kz*diag(dd_W_mean)*shear+Pr*Laplacian_square);
                                A12=O;
                                A13=Pr*Ra_T_q*Laplacian_inv*(-1i*kz*D1x);

                                A21=-1i*ky*diag(d_W_mean)*shear;
                                A22=-1i*kz*diag(W_mean)*shear+Pr*Laplacian;
                                A23=Pr*Ra_T_q*1i*ky*I;

                                A31=-diag(d_T_mean)+(1-flux_T*wT_int)*1i*kz*D1x/(ky^2+kz^2);
                                A32=(1-flux_T*wT_int)*(-1i*ky)/(ky^2+kz^2)*I;
                                A33=-1i*kz*diag(W_mean)*shear+Laplacian;

                                A=[A11,A12,A13;
                                    A21,A22,A23;
                                    A31,A32,A33];

                                [eig_vec,eig_val]=eig(A,M);
                                eig_val=diag(eig_val);
                                eig_val(find(real(eig_val)==Inf))=-Inf;
                                [growth_rate_no_shear_com(Lx_ind,n_elevator_ind,Ra_T_q_ind,Pr_ind,ky_ind,kz_ind),max_ind]=max(real(eig_val));

                            end

                        end
                        [~,kz_max_ind]=max(squeeze(growth_rate(Lx_ind,n_elevator_ind,Ra_T_q_ind,Pr_ind,ky_ind,:)));
                        kz_max(Lx_ind,n_elevator_ind,Ra_T_q_ind,Pr_ind,ky_ind)=kz_list(kz_max_ind);
                        [~,kz_0_ind]=min(squeeze(abs(growth_rate(Lx_ind,n_elevator_ind,Ra_T_q_ind,Pr_ind,ky_ind,3:end))));%This will get rid of the first two points that always gives growth rate close to zero 
                        kz_0(Lx_ind,n_elevator_ind,Ra_T_q_ind,Pr_ind,ky_ind)=kz_list(kz_0_ind);
                    end
                end
            end
        end
    end
    save(['flux_T_RBC_',flag.name,'.mat']);
else
    load(['flux_T_RBC_',flag.name,'.mat']);
end

switch flag.name
    case 'growth_rate'
        data{1}.x=kz_list;
        data{1}.y=squeeze(growth_rate);
        %data{2}.x=kz_list;
        %data{2}.y=squeeze(growth_rate_no_flux_com);
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k-','b--'};
        plot_config.linewidth=3;
        plot_config.label_list={1,'$k_z$','max[Re($\lambda$)]'};
        plot_config.print_size=[1,1000,900];
        plot_config.legend_list={0};
        %plot_config.legend_list={1,'With flux feedback','No flux feedback'};
        plot_config.fontsize_legend=28;
        plot_config.name=['flux_RBC_',flag.name,'_Ra_Tq=',num2str(Ra_T_q_list),...
            '_Pr=',num2str(Pr_list),'_Lx=',num2str(Lx_list),'_kx_mean=',num2str(kx_mean),...
            '_flux_T=',num2str(flux_T),'.png'];
        plot_line(data,plot_config);
    case 'Lx'
        data{1}.x=Lx_list(2:end);
        data{1}.y=squeeze(kz_0(2:end));
        data{2}.x=Lx_list(2:end);
        data{2}.y=squeeze(kz_max(2:end));
        data{3}.x=linspace(0,max(data{1}.x),100);
        data{3}.y=2*pi*ones(size(data{3}.x));
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k*','bo','r--'};
        plot_config.linewidth=2;
        plot_config.label_list={1,'$L_x$',''};
        plot_config.print_size=[1,1000,900];
        plot_config.legend_list={1,'$k_{z,0}$','$k_{z,max}$'};
        plot_config.fontsize_legend=36;
        plot_config.name=['flux_RBC_',flag.name,'_Ra_Tq=',num2str(Ra_T_q_list),...
            '_Pr=',num2str(Pr_list),...
            '_flux_T=',num2str(flux_T),'.png'];
        plot_line(data,plot_config);
    case 'Ra_T_q'
        data{1}.x=Ra_T_q_list(2:end);
        data{1}.y=squeeze(kz_0(2:end));
        data{2}.x=Ra_T_q_list(2:end);
        data{2}.y=squeeze(kz_max(2:end));
        data{3}.x=linspace(10000,max(data{1}.x),100);
        data{3}.y=2*pi*ones(size(data{3}.x));
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k*','bo','r--'};
        plot_config.linewidth=2;
        plot_config.loglog=[1,0];
        plot_config.label_list={1,'$Ra_{T,q}$',''};
        plot_config.print_size=[1,1000,900];
        plot_config.legend_list={1,'$k_{z,0}$','$k_{z,max}$'};
        plot_config.fontsize_legend=36;
        plot_config.ylim_list=[1,3,10];
        plot_config.name=['flux_RBC_',flag.name,...
            '_Pr=',num2str(Pr_list),'_Lx=',num2str(Lx),...
            '_flux_T=',num2str(flux_T),'.png'];
        plot_config.xtick_list=[1,10^4,10^5,10^6,10^7,10^8];
        plot_line(data,plot_config);
        
    case 'Pr'
        data{1}.x=Pr_list;
        data{1}.y=squeeze(kz_0);
        data{2}.x=Pr_list;
        data{2}.y=squeeze(kz_max);
        data{3}.x=linspace(min(data{1}.x),max(data{1}.x),100);
        data{3}.y=2*pi*ones(size(data{3}.x));
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k*','bo','r--'};
        plot_config.linewidth=2;
        plot_config.loglog=[1,0];
        plot_config.label_list={1,'$Pr$',''};
        plot_config.print_size=[1,1000,900];
        plot_config.legend_list={1,'$k_{z,0}$','$k_{z,max}$'};
        plot_config.fontsize_legend=36;
        plot_config.name=['flux_RBC_',flag.name,...
            '_Ra_Tq=',num2str(Ra_T_q_list),'_Lx=',num2str(Lx),...
            '_flux_T=',num2str(flux_T),'.png'];
        plot_config.xtick_list=[1,0.01,0.1,1,7,70,700];
        plot_line(data,plot_config);
    case 'ke'
        data{1}.x=n_elevator_list*2*pi/Lx_list(1);
        data{1}.y=squeeze(kz_0);
        data{2}.x=n_elevator_list*2*pi/Lx_list(1);
        data{2}.y=squeeze(kz_max);
        data{3}.x=linspace(0,max(data{1}.x),100);
        data{3}.y=2*pi*ones(size(data{3}.x));
        
        fit_kz_0=lsqr([data{1}.x',ones(6,1)],data{1}.y');
        data{4}.x=data{1}.x;
        data{4}.y=fit_kz_0(1)*data{1}.x+fit_kz_0(2);
        
        fit_kz_max=lsqr([data{2}.x',ones(6,1)],data{2}.y');
        data{5}.x=data{2}.x;
        data{5}.y=fit_kz_max(1)*data{5}.x+fit_kz_max(2);
        plot_config.Markerindex=3;
        plot_config.user_color_style_marker_list={'k*','bo','r--','k-.','b-.'};
        plot_config.linewidth=2;
        plot_config.label_list={1,'$k_e$',''};
        plot_config.print_size=[1,1000,900];
        plot_config.legend_list={1,'$k_{z,0}$','$k_{z,max}$'};
        plot_config.fontsize_legend=36;
        plot_config.name=['flux_RBC_',flag.name,'_Ra_Tq=',num2str(Ra_T_q_list),...
            '_Pr=',num2str(Pr_list),'_Lx=',num2str(Lx_list),...
            '_flux_T=',num2str(flux_T),'.png'];
        plot_line(data,plot_config);
end


error('1');

data{1}.x=kz_list_num;
data{1}.y=growth_rate(1,1,1,1,length(kz_list_int)+1:end);
% data{2}.x=kz_list_int;
% data{2}.y=growth_rate(1,1,1,1,1:length(kz_list_int));
plot_config.Markerindex=3;
plot_config.user_color_style_marker_list={'k-','bo'};
plot_config.label_list={1,'$k_z$','$max[Re(\lambda)]$'};
plot_config.print_size=[1,1000,900];
plot_config.name=['flux_RBC_growth_rate_Ra_Tq=',num2str(Ra_T_q),...
    '_Pr=',num2str(Pr),'_Lx=',num2str(Lx),'_kx_mean=',num2str(kx_mean),...
    '_flux_T=',num2str(flux_T),'.png'];
plot_line(data,plot_config);

