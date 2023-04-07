classdef settling_post
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
        Prandtl;
        tau;
        R_rho;
        W_st;
        A_noise;
        d_T_z;
        d_C_z;
        
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
        
        u_coeff;
        w_coeff;
        spectrum_TKE;
        
        title_time=1;
        
        post_store_dt;
        
    end
    
    methods
        function obj = settling_post(h5_name,flag)
            %dedalus_post Construct an instance of this class
            %   Detailed explanation goes here
            if nargin<2 || isempty(flag)
                flag.print=1;
                flag.video=1;
                flag.visible=1;
            end
            %construction function...
            obj.h5_name=h5_name;%name for the h5file
            obj.h5_name_scalar_data=strrep(h5_name,'analysis','scalar_data');
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
            obj.kz_list=h5read_complex(h5_name,'/scales/kz');        

        end
               
        
        
        function obj=snapshot(obj,variable_name,video_ind,zlim_list,snap_ind)
            %%plot the snapshot of salinity and generate video if any
            if nargin<3 || isempty(video_ind)
               video_ind=1; 
            end
                        
            obj.(variable_name)=h5read_complex(obj.h5_name,['/tasks/',variable_name]);

            plot_config.colormap='bluewhitered';%bluewhitered
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

                    plot_config.fontsize=40;
                    plot_config.ylim_list=[1,round(min(data{1}.y)),round(max(data{1}.y))];
                    %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
                    
                    if obj.title_time
                        plot_config.title_list={1,['$t=$',num2str(round(obj.t_list(t_ind)))]};
                    else
                        plot_config.title_list={0};
                    end
                    
                    plot_config.print_size=[1,900,1200];
                    plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_',num2str(round(obj.t_list(t_ind))),'.png'];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
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
                    %plot_config.ytick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0];

                    if obj.title_time
                        plot_config.title_list={1,['$t=$',num2str(round(obj.t_list(t_ind)))]};
                    else
                        plot_config.title_list={0};
                    end

                    plot_config.print_size=[1,900,1200];
                    plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_ind_',num2str(t_ind),'.png'];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
                    plot_config.fontsize=35;
                    plot_contour(data,plot_config);
                    %plot_config.label_list={1,'$x$',''};
                    %plot_config.name=[obj.h5_name(1:end-3),'_snapshot_',variable_name,'_t_',num2str(t_ind),'_no_ylabel.png'];
                    %plot_contour(data,plot_config);
                end
            end
        end
        
        
        
        function obj=spectrum_snapshot(obj,variable_name)
            %%plot the spectrum of salnity as time varies, also generate
            %%video if any
            obj.([variable_name,'_coeff'])=h5read_complex(obj.h5_name,['/tasks/',variable_name,'_coeff']);
            %coeff.r+1i*coeff.i;
            if obj.video
                for t_ind=1:length(obj.t_list)
                    clear data plot_config;
                    
                    data{1}.x=obj.kx_list;
                    data{1}.y=obj.kz_list(1:obj.Nz/2);
                    plot_config.label_list={1,'$k_x$','$k_z$'};
                    data{1}.z=log10(abs(obj.([variable_name,'_coeff'])(1:obj.Nz/2,:,t_ind)).^2);
                    %data{2}.x=obj.kx_list/obj.k_opt;
                    %data{2}.y=obj.ks/obj.k_opt*ones(size(obj.kx_list));
                    plot_config.print_size=[1,1100,900];
                    plot_config.loglog=[0,0];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
                    %plot_config.xtick_list=[0.01,0.1,1,10];

                    plot_config.colormap='white_zero';
                    plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_t_',num2str(round(obj.t_list(t_ind),2)),'.png'];
                    frame_spectrum_2D(t_ind)=plot_contour(data,plot_config);

                    dx=diff(obj.kx_list); dx=dx(1);
                    dz=diff(obj.kz_list); dz=dz(1);
                    
                    data{1}.x=obj.kx_list;
                    data{2}.x=obj.kz_list(1:obj.Nz/2);
                    plot_config.label_list={1,'$k_x$ or $k_z$',''};
                    data{1}.y=2*dz*sum(abs(obj.([variable_name,'_coeff'])(1:obj.Nz/2,:,t_ind)).^2,1);
                    data{2}.y=2*dx*sum(abs(obj.([variable_name,'_coeff'])(1:obj.Nz/2,:,t_ind)).^2,2);
                    
                    plot_config.loglog=[1,1];
                    plot_config.ytick_list=[0,0.001,0.01,0.1,1,10,100,1000];
                    plot_config.ylim_list=[0];%,0.1,10];
                    plot_config.xtick_list=[1,0.001,0.01,0.1,1,10,100];
                    switch variable_name
                        case 'forcing_var_x'
                            variable_name_legend='{f_x}';
                        case 'forcing_var_y'
                            variable_name_legend='{f_y}';
                        otherwise
                            variable_name_legend=variable_name;
                    end
                    plot_config.legend_list={1,['$\int E_',variable_name_legend,'(k_x,k_z)dk_z$'],['$\int E_',variable_name_legend,'(k_x,k_z)d k_x$']};
                    plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_1D_t_',num2str(round(obj.t_list(t_ind),2)),'.png'];
                    plot_config.print=obj.print;
                    plot_config.visible=obj.visible;
                    frame_spectrum_1D(t_ind)=plot_line(data,plot_config);
                end
               plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_2D_t_video.mp4'];
               plot_video(frame_spectrum_2D,plot_config);
               plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_1D_t_video.mp4'];
               plot_video(frame_spectrum_1D,plot_config);
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
            switch variable_name
                case 'forcing_var_x'
                    variable_name_legend='{f_x}';
                case 'forcing_var_y'
                    variable_name_legend='{f_y}';
                otherwise
                    variable_name_legend=variable_name;
            end
            plot_config.legend_list={1,['$\int E_',variable_name_legend,'(k_x,k_z)dk_z$'],['$\int E_',variable_name_legend,'(k_x,k_z)d k_x$']};
            plot_config.name=[obj.h5_name(1:end-3),'_spectrum_',variable_name,'_1D_time_average.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_config.ytick_list=[0,10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),...
                0.001,0.01,0.1,1,10,100,1000];
            plot_config.xtick_list=[1,10^(-8),10^(-7),10^(-6),10^(-5),10^(-4),...
                0.001,0.01,0.1,1,10,100,1000];
            plot_line(data,plot_config);
        end
        
        
        function obj=spectrum_TKE_average(obj)
            %%This function plot the kinetic energy, (deduct the background laminar flow)
            %%This is average over time... plot the turbulence kinetic
            %%energy spectrum as kx and kz
            
            %%This is the post-processing for the TKE in 2D... 
            obj.w_coeff=h5read_complex(obj.h5_name,'/tasks/w_coeff');
            %=w_coeff.r+1i*w_coeff.i;
            obj.u_coeff=h5read_complex(obj.h5_name,'/tasks/u_coeff');
            %obj.u_coeff=u_coeff.r+1i*u_coeff.i;
            for t_ind=1:length(obj.t_list)
                obj.spectrum_TKE(:,:,t_ind)=abs(obj.u_coeff(:,:,t_ind)).^2+abs(obj.w_coeff(:,:,t_ind)).^2;
            end
            spectrum_TKE_average=mean(abs(obj.spectrum_TKE(1:obj.Nz/2,:,1:end)).^2,3);%max(3*max_ind,30)

            data{1}.z=log10(spectrum_TKE_average);
            data{1}.x=obj.kx_list;
            data{1}.y=obj.kz_list(1:obj.Nz/2);
            
            plot_config.loglog=[0,0];
            plot_config.print_size=[1,1100,900];
            plot_config.label_list={1,'$k_x$','$k_z$'};
            plot_config.colormap='white_zero';
            plot_config.name=[obj.h5_name(1:end-3),'_spectrum_TKE_2D_time_average.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_contour(data,plot_config);
            
            dx=diff(obj.kx_list); dx=dx(1);
            dz=diff(obj.kz_list); dz=dz(1);

            data{1}.x=obj.kx_list;
            data{1}.y=2*dz*sum(spectrum_TKE_average,1);
            data{2}.x=obj.kz_list(1:obj.Nz/2);
            data{2}.y=2*dx*sum(spectrum_TKE_average,2);
            
            plot_config.loglog=[1,1];
            plot_config.ytick_list=[1,0.001,0.01,0.1,1,10,100,1000];
            plot_config.label_list={1,'$k_x$ or $k_z$',''};
            plot_config.legend_list={1,'$\int E_u(k_x,k_z)dk_z$','$\int E_u(k_x,k_z)d k_x$'};
            plot_config.name=[obj.h5_name(1:end-3),'_spectrum_TKE_1D_time_average.png'];
            plot_config.print=obj.print;
            plot_config.visible=obj.visible;
            plot_line(data,plot_config);
        end
        
        function obj=scalar_data(obj,variable_name)
            obj.([variable_name])=h5read_complex(obj.h5_name_scalar_data,['/tasks/',variable_name]);
            data{1}.x=obj.t_list;
            data{1}.y=squeeze(obj.([variable_name]));
            plot_config.label_list={1,'$t$',variable_name};
            plot_config.name=[obj.h5_name(1:end-3),variable_name,'_time_history.png'];
            plot_line(data,plot_config);
        end
        
    end
end

