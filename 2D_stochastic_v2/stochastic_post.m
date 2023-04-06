classdef stochastic_post
    %POST_CLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x_list;
        y_list;
        kx_list;
        ky_list;
        t_list;
        
        Nx;
        Ny;
        eps;
        k1;
        k2;
        nu;
        stop_sim_time;
        initial_dt;
        Lx;
        Ly;
       
        h5_name;
        print;
        visible;
        no_ylabel;
        video;
        
    end
    
    methods
        function obj = stochastic_post(h5_name,flag)
            %dedalus_post Construct an instance of this class
            %   Detailed explanation goes here
            if nargin<2 || isempty(flag)
                flag.print=1;
                flag.video=1;
                flag.visible=1;
            end
            %construction function...
            obj.h5_name=h5_name;%name for the h5file
            
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
            obj.y_list=h5read(h5_name,'/scales/y/1.0');
            obj.t_list=h5read_complex(h5_name,'/scales/sim_time');
            obj.kx_list=h5read_complex(h5_name,'/scales/kx');        
            obj.ky_list=h5read_complex(h5_name,'/scales/ky');        

        end
               
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

