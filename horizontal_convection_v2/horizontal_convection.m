classdef horizontal_convection_post
    %POST_CLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x_list;
        z_list;
        kx_list;
        kz_list;
        t_list;
        
        Nx;
        Nz;
        eps;
        Rayleigh;
        Prandtl;
        tau;
        R_rho;
        W_st;
        A_noise;
        d_T_z;
        d_C_z;
        
        modulation;
        
        stop_sim_time;
        initial_dt;
        Lx;
        Lz;
       
        h5_name;
        print;
        visible;
        no_ylabel;
        video;
        
        u;
        w;
        T;
        T_tot;

        u_coeff;
        w_coeff;
        
        spectrum_TKE;
        
        title_time=1;
        
        post_store_dt;
        x_Test;
        
        z_basis_mode='Chebyshev';
        spec_kx_z_b;
        
        phi;
        kappa;
        A_LS;
        gaussian_sigma;
        restart_t0;
        
        spec_kx_z_T;
    end
    
    methods
        function obj = horizontal_convection_post(h5_name,flag)
            %dedalus_post Construct an instance of this class
            %   Detailed explanation goes here
            if nargin<2 || isempty(flag)
                flag.print=1;
                flag.video=1;
                flag.visible=1;
            end
            %construction function...
            obj.h5_name=h5_name;%name for the h5file
            %obj.h5_name_scalar_data=strrep(h5_name,'analysis','scalar_data');
            %modify these flag.
            obj.print=flag.print;
            obj.video=flag.video;
            obj.visible=flag.visible;
            obj.no_ylabel=flag.no_ylabel;
            %display h5 file
            %h5disp(h5_name);
            
            %read the flag_table.
            flag_table=readtable([h5_name(1:end-14),'flag.txt']);
            for table_ind=1:length(flag_table.x_Test)
               if isnumeric(obj.(flag_table.x_Test{table_ind}(3:end)))
                   obj.(flag_table.x_Test{table_ind}(3:end))=str2num(flag_table.x123_{table_ind}(1:end-1));
               else
                   obj.(flag_table.x_Test{table_ind}(3:end))=flag_table.x123_{table_ind}(1:end-1);
               end
            end
            
            %POST_CLASS Construct an instance of this class
            %   Detailed explanation goes here
            obj.x_list=h5read(h5_name,'/scales/x/1.0');
            obj.z_list=h5read(h5_name,'/scales/z/1.0');
            obj.t_list=h5read_complex(h5_name,'/scales/sim_time');
            obj.kx_list=h5read_complex(h5_name,'/scales/kx');        

        end
               
        
        
        function obj=snapshot(obj,variable_name,video_ind,zlim_list,snap_ind)
            %%plot the snapshot of salinity and generate video if any
            if nargin<3 || isempty(video_ind)
               video_ind=1; 
            end
            if strcmp(variable_name,'T_tot')            
                T=h5read_complex(obj.h5_name,['/tasks/',variable_name(1)]);
                
                background_T=(obj.kappa-1)*obj.z_list*ones(size(obj.x_list))'+1;
                for t_ind=1:length(obj.t_list)
                    T_tot(:,:,t_ind)=T(:,:,t_ind)+background_T;
                end
                obj.(variable_name)=T_tot;
                plot_config.colormap='jet';
            else
                obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                plot_config.colormap='bluewhitered';%bluewhitered
    
            end
            variable_max=max(max(max(obj.(variable_name))));
            variable_min=min(min(min(obj.(variable_name))));
            plot_config.zlim_list=[1,variable_min,variable_max];

            
            %Update 2022/09/14, update the zlim_list option so I could turn
            %off this. 
            if nargin<4 || isempty(zlim_list)
                %do nothing
            else
               plot_config.zlim_list=zlim_list; 
            end
            
            if obj.video
                snapshot_ind=1;
                for t_ind=1:video_ind:length(obj.t_list)
                    data{1}.z=obj.(variable_name)(:,:,t_ind);
                    
                    data{1}.x=obj.x_list;
                    data{1}.y=obj.z_list;
                    plot_config.label_list={1,'$x$','$z$'};

                    plot_config.fontsize=28;
                    %plot_config.ylim_list=[1,round(min(data{1}.y)),round(max(data{1}.y))];
                    %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
                    
                    if obj.title_time
                        plot_config.title_list={1,['$t=$',num2str(round(obj.t_list(t_ind)))]};
                    else
                        plot_config.title_list={0};
                    end
                    
                    plot_config.print_size=[1,2000,500];
                    plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_',num2str(round(obj.t_list(t_ind))),'.png'];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
                    plot_config.ylim_list=[1,round(min(data{1}.y)),round(max(data{1}.y))];
                    plot_config.xlim_list=[1,round(min(data{1}.x)),round(max(data{1}.x))];
                    plot_config.zlim_list=[1,round(min(min(data{1}.z))),round(max(max(data{1}.z)))];
                    plot_config.ztick_list=[1,0,0.5,1];
                    plot_config.daspect=[1,1,1,1];
                    plot_config.fontsize=28;
                    snapshot(snapshot_ind)=plot_contour(data,plot_config);
                    snapshot_ind=snapshot_ind+1;
                    %plot_config.label_list={1,'$x$',''};
                    %plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_',num2str(round(obj.t_list(t_ind))),'_no_ylabel.png'];
                    %plot_contour(data,plot_config);
                end
               plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_video.mp4'];
               plot_video(snapshot,plot_config);
            end
            
            if nargin<5 || isempty(snap_ind)
                %do nothing
            else
                for t_ind=snap_ind
                    data{1}.z=obj.(variable_name)(:,:,t_ind);

                    data{1}.x=obj.x_list;
                    data{1}.y=obj.z_list;
                    plot_config.label_list={1,'$x$','$z$'};
            
%                     plot_config.fontsize=40;
                    plot_config.ylim_list=[1,round(min(data{1}.y)),round(max(data{1}.y))];
                    plot_config.xlim_list=[1,round(min(data{1}.x)),round(max(data{1}.x))];
                    plot_config.zlim_list=[1,round(min(min(data{1}.z))),round(max(max(data{1}.z)))];
                    plot_config.ztick_list=[1,0,0.5,1];
                    %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0];

                    if obj.title_time
                        plot_config.title_list={1,['$t=$',num2str(round(obj.t_list(t_ind)))]};
                    else
                        plot_config.title_list={0};
                    end
                    plot_config.print_size=[1,2000,500];
                    plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_ind_',num2str(t_ind),'.png'];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
                    plot_config.fontsize=28;
                    plot_config.daspect=[1,1,1,1];
                    plot_contour(data,plot_config);
                    %plot_config.label_list={1,'$x$',''};
                    %plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_',num2str(t_ind),'_no_ylabel.png'];
                    %plot_contour(data,plot_config);
                end
            end
        end
        
        
        
        function obj=spectrum_average(obj,variable_name,t_range)
            %%This function plot the 
            %%plot the overall spectrum averaged over time
            if nargin<3 || isempty(t_range)
               t_range(1)=obj.t_list(1);
               t_range(2)=obj.t_list(end); 
            end
           
            [val,t_ind_begin]=min(abs(obj.t_list-t_range(1)));
            [val,t_ind_end]=min(abs(obj.t_list-t_range(2)));
            
            %coeff.r+1i*coeff.i;
            
%             obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
%             for t_ind=1:length(obj.t_list)
%                 obj.(['E_',variable_name])(t_ind)=sum(sum(obj.(variable_name)(:,:,t_ind).^2))/obj.Nx/obj.Nz/2;
%             end
%             [val,max_ind]=max(obj.(['E_',variable_name]));
%            
            if strcmp(obj.z_basis_mode,'Chebyshev')
                %modify for the vertical as the bounded domain.. Chebyshev
                %basis instead of the Fourier.. so we do not have spectrum
                %int he vertical
                obj.([variable_name])=h5read_complex(obj.h5_name,['/tasks/',variable_name]);

                x_list=obj.x_list;
                dx=mean(diff(obj.x_list));
                Nx=length(x_list);
                Fs=1/dx;
                obj.kx_list=2*pi*Fs*(0:(Nx/2-1))/Nx; %Note that this wavenumber corresponding to circular frequency
                t_ave_num=0;
                spec_kx_z=zeros(length(obj.kx_list),length(obj.z_list));
                for t_ind=t_ind_begin:t_ind_end
                    for z_ind=1:length(obj.z_list)
                        spec_tmp=abs(fft(obj.([variable_name])(z_ind,:,t_ind))/Nx);
                        spec_tmp=spec_tmp(1:Nx/2);
                        %spec_tmp(2:end)=2*spec_tmp(2:end);
                        
                        spec_kx_z(:,z_ind)=spec_kx_z(:,z_ind)+(spec_tmp)';
                    end
                    t_ave_num=t_ave_num+1;
                end
                
                spec_kx_z=spec_kx_z/t_ave_num;
%                 spec_kx_z=abs(spec_kx_z);
                obj.(['spec_kx_z_',variable_name])=spec_kx_z;
                data{1}.x=obj.kx_list;
                data{1}.y=obj.z_list;
                data{1}.z=(spec_kx_z)';
                plot_config.label_list={1,'$k_x$','$z$'};

                %spectrum_average=mean(abs(obj.([variable_name,'_coeff'])(:,:,t_ind_begin:t_ind_end)).^2,3);
                %data{1}.z=log10(spectrum_average);
                plot_config.loglog=[0,0];
                plot_config.print_size=[1,800,900];
                plot_config.colormap='white_zero';
                plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_time_average.png'];
                plot_config.print=obj.print;
                plot_config.ytick_list=[0,0,0.2,0.4,0.6,0.8,1];
                plot_config.xtick_list=[1,0,10,20,30];
                plot_config.xlim_list=[1,0,20];
                %plot_config.ylim_list=[1,0,1];
                
                %This range is for the plotting of the bounded salt-finger
                plot_config.zlim_list=[1,0,0.16];
                %plot_config.ztick_list=[1,0,0.04,0.08,0.12,0.16]
                
                plot_config.zlim_list=0;
                plot_config.ztick_list=0;
                plot_config.xlim_list=[1,0,50];
                
                plot_config.visible=obj.visible;
                plot_config.fontsize=32;
                plot_contour(data,plot_config);

            elseif strcmp(obj.z_basis_mode,'Fourier')
                obj.([variable_name,'_coeff'])=h5read_complex(obj.h5_name,['/tasks/',variable_name,'_coeff']);

                data{1}.x=obj.kx_list;
                data{1}.y=obj.kz_list(1:obj.Nz/2);
                plot_config.label_list={1,'$k_x$','$k_z$'};

                spectrum_average=mean(abs(obj.([variable_name,'_coeff'])(1:obj.Nz/2,:,1:end)).^2,3);
                data{1}.z=log10(spectrum_average);
                plot_config.loglog=[0,0];
                plot_config.print_size=[1,1100,900];
                plot_config.colormap='white_zero';
                plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_time_average.png'];
                plot_config.print=obj.print;
                plot_config.ytick_list=[0,10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),...
                    0.001,0.01,0.1,1,10,100,1000];
                plot_config.xtick_list=[0,10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),...
                    0.001,0.01,0.1,1,10,100,1000];
                plot_config.visible=obj.visible;
                plot_contour(data,plot_config);

                dx=diff(obj.kx_list); dx=dx(1);
                dz=diff(obj.kz_list); dz=dz(1);

                data{1}.x=obj.kx_list;
                data{1}.y=2*dz*sum(spectrum_average,1);
                data{2}.x=obj.kz_list(1:obj.Nz/2);
                data{2}.y=2*dx*sum(spectrum_average,2);

                plot_config.loglog=[1,1];

                plot_config.label_list={1,'$k_x$ or $k_z$',''};
                plot_config.legend_list={1,['$\int E_',variable_name,'(k_x,k_z)dk_z$'],['$\int E_',variable_name,'(k_x,k_z)d k_x$']};
                plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_1D_time_average.png'];
                plot_config.print=obj.print;
                plot_config.visible=obj.visible;
                plot_line(data,plot_config);
            end
        end
        
        function obj=x_ave(obj,variable_name,t_ind)
            %%plot the streamwise averaged salnity
            %%as a function of z (vertical axis) and time
            if nargin<3 || isempty(t_ind)
                %The default option, just average over the second half of
                %data...
                t_ind_begin=1;
                t_ind_end=length(obj.t_list);
            else
                t_ind_begin=t_ind(1);
                if length(t_ind)==1
                    t_ind_end=length(obj.t_list);
                else
                    t_ind_end=t_ind(2);
                end
            end
            
            data{1}.x=obj.t_list(t_ind_begin:t_ind_end);
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'$t$','$z/l_{opt}$'};
%             else
            data{1}.y=obj.z_list;
            plot_config.label_list={1,'$t$','$z$'};
%             end
            if strcmp(obj.flow,'HB_benard_shear_periodic')
                switch variable_name
                    case {'u','w','T'}
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',upper(variable_name),'_0']);
                        data{1}.z=real(obj.(variable_name)(:,t_ind_begin:t_ind_end));
                    case 'dy_T_mean_q'
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                        data{1}.z=obj.(variable_name)(:,t_ind_begin:t_ind_end);
                end
            else
                switch variable_name
                    case {'u','v','w','S','T','p','dy_T_mean_q','dy_S_mean_q'}
                        obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                    case {'uS','wS','uT','wT','uw','ww'}%%
                        var_1=variable_name(1);
                        var_2=variable_name(2);
                        if strcmp(var_1,'u') %%this require in default, the u is always in the first variable....
                           obj=obj.u_fluctuation_read();
                           var_1_data=obj.u_fluctuation;
                        else
                           var_1_data=h5read_complex(obj.h5_name,['/tasks/',var_1]);
                        end
                        var_2_data=h5read_complex(obj.h5_name,['/tasks/',var_2]);
                        obj.(variable_name)=var_1_data.*var_2_data;
                end

                data{1}.z=squeeze(mean(obj.(variable_name)(:,:,t_ind_begin:t_ind_end),2));

            end
            %             plot_config.label_list={1,'$t$','$z/l_{opt}$'};
            plot_config.colormap='bluewhitered';
            plot_config.print_size=[1,1600,1600];
%             plot_config.ztick_list=[1,-0.001,0.001];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_x_ave.png'];
            plot_config.ylim_list=[1,round(min(data{1}.y),1),round(max(data{1}.y),1)];
            if round(min(data{1}.x),1)==round(max(data{1}.x),1)
                plot_config.xlim_list=[1,min(data{1}.x),max(data{1}.x)];
            else
                plot_config.xlim_list=[1,round(min(data{1}.x),1),round(max(data{1}.x),1)];
            end
            plot_config.ytick_list=[1,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
            plot_config.fontsize=28;
            plot_contour(data,plot_config);
            
            if strcmp(variable_name,'u')
                for t_ind=1:length(data{1}.x)
                    [val,ind]=mink(abs(data{1}.z(:,t_ind)),5);
                    
                    Nz=obj.Nz;
                    if abs(mod(ind(2),Nz)-mod(ind(1),Nz))>1
                        z_max(t_ind)=(data{1}.y(ind(1))+data{1}.y(ind(2)))/2;
                    elseif abs(mod(ind(3),Nz)-mod(ind(1),Nz))>2
                        z_max(t_ind)=(data{1}.y(ind(1))+data{1}.y(ind(3)))/2;
                    else
                        z_max(t_ind)=(data{1}.y(ind(1))+data{1}.y(ind(4)))/2;
                    end
                end
                
                data_z{1}.x=data{1}.x;
                data_z{1}.y=z_max';
                plot_config.label_list={1,'$t$','$z$'};
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_TW_c_z.png'];
                plot_config.ylim_list=0;
                plot_line(data_z,plot_config);
                
                for t_ind=1:length(z_max)-5
                    if abs(z_max(t_ind+1)-z_max(t_ind))>0.2
                         if all(abs(z_max(t_ind+2:t_ind+5)-z_max(t_ind+1))<0.2)
                             z_max(t_ind+1:end)=z_max(t_ind+1:end)+0.5;
                         else
                             z_max(t_ind+1)=z_max(t_ind+1)+0.5;
                         end
%                         if abs(z_max(t_ind+2)-z_max(t_ind+1))>0.45
%                             z_max(t_ind+1)=z_max(t_ind+1)+0.5;
%                         else
%                             z_max(t_ind+1:end)=z_max(t_ind+1:end)+0.5;
%                         end
                    end
                end
                data_z{1}.y=z_max';
                plot_config.label_list={1,'$t$','$z$'};
                plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_TW_c_z_processed.png'];
                plot_config.ytick_list=0;
                plot_line(data_z,plot_config);
                
                coeff=[ones(length(data_z{1}.x),1),data_z{1}.x]\z_max';
                obj.phase_c_z=coeff(2);
                
            end
            
        end
        
        function obj=total_xt_ave(obj,variable_name,x_ind,t_ind)
%             data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%             dz=diff(obj.z_list); dz=dz(1);
%             z_list_full=[obj.z_list;obj.z_list(end)+dz];
            x_len=length(obj.x_list);
            time_len=length(obj.t_list);
            
            if nargin<4 || isempty(x_ind)
                %The default option, just average over the second half of
                %data...
                x_ind_begin=1;
                x_ind_end=x_len;
            else
                x_ind_begin=x_ind(1);
                if length(x_ind)==1
                    x_ind_end=x_len;
                else
                    x_ind_end=x_ind(2);
                end
            end
            
            
            if nargin<4 || isempty(t_ind)
                %The default option, just average over the second half of
                %data...
                t_ind_begin=time_len/2;
                t_ind_end=time_len;
            else
                t_ind_begin=t_ind(1);
                if length(t_ind)==1
                    t_ind_end=time_len;
                else
                    t_ind_end=t_ind(2);
                end
            end

            switch variable_name
                case {'T','S'}
                    if obj.(['flux_',variable_name])
                        if strcmp(obj.flow,'HB_benard_shear_periodic')
                            variable_data=h5read_complex(obj.h5_name,['/tasks/',variable_name,'_0']);
                            data{2}.x=mean(obj.(['dy_',variable_name,'_mean_q'])(1,t_ind_begin:t_ind_end),2)*obj.z_list;
                            data{1}.x=obj.(['Pe_',variable_name])*squeeze(mean(variable_data(:,t_ind_begin:t_ind_end),2))+mean(obj.(['dy_',variable_name,'_mean_q'])(1,t_ind_begin:t_ind_end),2)*obj.z_list;
                        else
                            variable_data=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                            data{2}.x=mean(mean(obj.(['dy_',variable_name,'_mean_q'])(1,1,t_ind_begin:t_ind_end),2),3)*obj.z_list;
                            data{1}.x=obj.(['Pe_',variable_name])*squeeze(mean(mean(variable_data(:,:,t_ind_begin:t_ind_end),2),3))+mean(mean(obj.(['dy_',variable_name,'_mean_q'])(1,1,t_ind_begin:t_ind_end),2),3)*obj.z_list;
                        end
                    else
                        variable_data=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                        data{2}.x=obj.(['dy_',variable_name,'_mean'])*obj.z_list;
                        data{1}.x=obj.(['Pe_',variable_name])*squeeze(mean(mean(variable_data(:,x_ind_begin:x_ind_end,t_ind_begin:t_ind_end),2),3))+obj.(['dy_',variable_name,'_mean'])*obj.z_list;
                    
                    end
                    
                    if obj.(['dy_',variable_name,'_mean'])>0
                        plot_config.legend_list={1,['$ z+\langle ',variable_name,'\rangle_{h,t}$'],['$z$']};
                    else obj.(['dy_',variable_name,'_mean'])<0
                        data{2}.x=data{2}.x+1;
                        data{1}.x=data{1}.x+1;
                        plot_config.legend_list={1,['$1-z+\langle ',variable_name,'\rangle_{h,t}$'],['$1-z$']};
                        if obj.(['flux_',variable_name])
                            plot_config.legend_list={1,['$1+\langle\bar{\mathcal{T}}_{z,q}\rangle_t z+\langle ',variable_name,'\rangle_{h,t}$'],['$1+\langle\bar{\mathcal{T}}_{z,q}\rangle_t z$']};
                        end
                    end
                case {'rho'}
                    %error('not ready');
                    variable_data_T=h5read_complex(obj.h5_name,['/tasks/T']);
                    variable_data_S=h5read_complex(obj.h5_name,['/tasks/S']);

                    R_rho_T2S=obj.Ra_T/obj.Ra_S2T;
                    data{1}.x=-obj.dy_T_mean*obj.z_list+1/R_rho_T2S*obj.dy_S_mean*obj.z_list;
                    data{2}.x=-(obj.Pe_T*squeeze(mean(mean(variable_data_T,2),3))+obj.dy_T_mean*obj.z_list)...
                        +1/R_rho_T2S*(obj.Pe_S*squeeze(mean(mean(variable_data_S,2),3))+obj.dy_S_mean*obj.z_list);
                    plot_config.legend_list={1,['$-(\bar{\mathcal{T}}_z z+Pe_T \langle ','T','\rangle_h)+R_\rho^{-1}(\bar{\mathcal{S}}_z z+Pe_S \langle ','S','\rangle_h)$'],['$-\bar{\mathcal{T}}_z z+R_\rho^{-1}\bar{\mathcal{S}}_z z$']};
                    plot_config.fontsize_legend=24;
                case 'u'
                    u=h5read_complex(obj.h5_name,'/tasks/u');
                    syms z;
                    u_laminar=obj.F_sin/obj.ks^2*sin(obj.ks*z)...
                              +obj.F_sin_2ks/(2*obj.ks)^2*sin(2*obj.ks*z+obj.phase_2ks)...
                              +obj.F_sin_3ks/(3*obj.ks)^2*sin(3*obj.ks*z+obj.phase_3ks)...
                              +obj.F_sin_4ks/(4*obj.ks)^2*sin(4*obj.ks*z+obj.phase_4ks);
                    u_laminar_num=double(subs(u_laminar,z,obj.z_list));
                    data{1}.x=u_laminar_num;
                    data{2}.x=squeeze(mean(mean(u,2),3));
                    plot_config.legend_list={1,['$\bar{',variable_name,'}$'],['$\bar{',variable_name,'} +''\langle ',variable_name,'''\rangle_h$']};

            end
%             if strcmp(obj.flow(1:7),'IFSC_2D')
%                 data{1}.y=obj.z_list/(2*pi/obj.k_opt);
%                 data{2}.y=obj.z_list/(2*pi/obj.k_opt);
%                 plot_config.label_list={1,'','$z/l_{opt}$'};
%             else
                data{1}.y=obj.z_list;
                data{2}.y=obj.z_list;
%             end
            plot_config.label_list={1,'','$z$'};
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1];
            plot_config.ylim_list=[1,round(min(data{1}.y)),round(max(data{1}.y))];
%             plot_config.label_list={1,'','$z/l_{opt}$'};
            plot_config.print_size=[1,1600,1600];
            plot_config.print=obj.print;
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_xt_ave.png'];
            plot_config.linewidth=3;
            plot_line(data,plot_config);
            
            %data{1}.x=NaN; data{1}.y=NaN;
            %plot_config.label_list={1,plot_config.legend_list{3},'$z$'};
            %plot_config.legend_list={0};
            plot_config.print_size=[1,500,900];
            plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_xt_ave_profile_only.png'];
            plot_config.fontsize_legend=16;
            plot_config.linewidth=3;
            plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
            plot_line(data,plot_config);
            
            plot_config.print=obj.print;
            plot_config.visible=0;
            plot_config.print_size=[1,500,900];

            switch variable_name
                case {'T','S'}
                if obj.video

                    for t_ind=1:length(obj.t_list)
                        data{1}.x=squeeze(mean(variable_data(:,:,t_ind),2))+obj.z_list;
                        data{2}.x=obj.z_list;
                        plot_config.legend_list={1,['$ z+\langle ',variable_name,'\rangle_{h}$'],['$z$']};
                        if obj.title_time
                            plot_config.title_list={1,['$t=$',num2str(round(obj.t_list(t_ind)))]};
                        else
                            plot_config.title_list={0};
                        end
                        if obj.no_ylabel
                           plot_config.label_list={1,'',''}; 
                        end
                        plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_x_ave_',num2str(round(obj.t_list(t_ind))),'.png'];
                        snapshot(t_ind)=plot_line(data,plot_config);
                    end
                   plot_config.name=[obj.h5_name(1:end-3),'_',variable_name,'_total_xt_ave_video.mp4'];
                   plot_video(snapshot,plot_config);
                end
                otherwise
                    warning('No video for this case');
            end
            
        end
        
        function obj=get_Nu(obj,variable_name,t_ind)
            obj.t_list=h5read_complex(obj.h5_name,'/scales/sim_time');
            time_len=length(obj.t_list);
            if nargin<3 || isempty(t_ind)
                %The default option, just average over the second half of
                %data...
                t_ind_begin=time_len/2;
                t_ind_end=time_len;
            else
                t_ind_begin=t_ind(1);
                if length(t_ind)==1
                    t_ind_end=time_len;
                else
                    t_ind_end=t_ind(2);
                end
            end
            
            if obj.flux_T~=0 && obj.flux_S~=0
                if strcmp(obj.store_variable,'all')
                    d_variable_data=h5read_complex(obj.h5_name,['/tasks/d_',variable_name]);
                else
                    error('Not ready, this interpolation is wrong!!');
                    z_list_cheb=obj.z_list*2-1;
                    variable_data=h5read_complex(obj.h5_name,['/tasks/',variable_name]);
                    [x,DM]=chebdif(length(z_list_cheb),1);
                    D1=DM(:,:,1);
                    for x_ind=1:size(variable_data,2)
                        for t_ind=1:size(variable_data,3)
                            variable_data_y=squeeze(variable_data(:,x_ind,t_ind));
                            variable_data_y_int=chebint(variable_data_y,z_list_cheb);
                            d_variable_data(:,x_ind,t_ind)=D1*variable_data_y_int;
                        end
                    end
                end

                d_variable_data_total_xt_ave=obj.(['Pe_',variable_name])*squeeze(mean(mean(d_variable_data(:,:,t_ind_begin:t_ind_end),2),3))+obj.(['dy_',variable_name,'_mean']);
                obj.(['d_',variable_name])=d_variable_data;
                switch variable_name
                    case 'T'
                        obj.Nu=d_variable_data_total_xt_ave;
                    case 'S'
                        obj.Nu_S=d_variable_data_total_xt_ave;
                end
            
            else
                %Update 2023/04/07, take the average over the first and the
                %end of the local maximum. 
                if strcmp(obj.flow,'HB_benard_shear_periodic')
                    obj.Nu_T_t=(-1)./(squeeze(obj.(['dy_',variable_name,'_mean_q'])(1,t_ind_begin:t_ind_end)));
                    Nu_mid=(max(obj.Nu_T_t)+min(obj.Nu_T_t))/2;
                    ind_local_min=find(islocalmin(obj.Nu_T_t).*(obj.Nu_T_t<Nu_mid));
                    obj.Nu_T_t=obj.Nu_T_t(ind_local_min(1):ind_local_min(end));
                    %d_variable_data_total_xt_ave=mean((-1)./(squeeze(obj.(['dy_',variable_name,'_mean_q'])(:,t_ind_begin:t_ind_end))));
                    d_variable_data_total_xt_ave=mean(obj.Nu_T_t);
                else
                    obj.Nu_T_t=(-1)./(squeeze(obj.(['dy_',variable_name,'_mean_q'])(1,1,t_ind_begin:t_ind_end)));
                    Nu_mid=(max(obj.Nu_T_t)+min(obj.Nu_T_t))/2;
                    ind_local_min=find(islocalmin(obj.Nu_T_t).*(obj.Nu_T_t<Nu_mid));
                    obj.Nu_T_t=obj.Nu_T_t(ind_local_min(1):ind_local_min(end));
                    %d_variable_data_total_xt_ave=mean((-1)./(squeeze(obj.(['dy_',variable_name,'_mean_q'])(:,:,t_ind_begin:t_ind_end))));
                    d_variable_data_total_xt_ave=mean(obj.Nu_T_t);
                end
                 switch variable_name
                    case 'T'
                        obj.Nu=d_variable_data_total_xt_ave;
                        Ra_T_q=obj.Ra_T;
                        obj.Ra_T_no_q=Ra_T_q/mean(obj.Nu);
                    case 'S'
                        obj.Nu_S=d_variable_data_total_xt_ave;
                        Ra_S2T_q=obj.Ra_S2T;
                        obj.Ra_S2T_no_q=Ra_S2T_q/mean(obj.Nu_S);
                 end
                 
            end
            data{1}.x=obj.t_list(t_ind_begin:t_ind_end);
            
            %Update 2023/04/07, take from one local max to the end of local
            %maximum. 
            data{1}.x=data{1}.x(ind_local_min(1):ind_local_min(end));
            data{1}.y=obj.Nu_T_t;
            plot_config.label_list={1,'$t$','$nu(t)$'};
            plot_config.name=[obj.h5_name(1:end-3),'Nu_T_t.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_line(data,plot_config);
            
        end
        
        
    end
end

