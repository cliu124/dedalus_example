
clear all;
close all;
slurm_num={'20230610162211'};
Pr=1;
Ra=1e5;
P=inv(sqrt(Pr*Ra));
obj.h5_name=['D:\Data\dedalus_example\fixed_flux_RBC_v2\analysis_',...
        slurm_num{1},'\analysis_s1.h5'];
obj.b=h5read_complex(obj.h5_name,['/tasks/','b']);
obj.bz=h5read_complex(obj.h5_name,['/tasks/','bz']);
obj.w=h5read_complex(obj.h5_name,['/tasks/','w']);
obj.z_list=h5read(obj.h5_name,'/scales/z/1.0');
flux=squeeze(mean(obj.w.*obj.b-P*obj.bz,2));

T_mean=squeeze(mean(mean(obj.b,1),2));
gradient=squeeze(mean(obj.b(end,:,:),2)-mean(obj.b(1,:,:),2));
for t_ind=1:size(obj.b,3)
    b_fluct(:,:,t_ind)=obj.b(:,:,t_ind)-T_mean(t_ind)-gradient(t_ind)*obj.z_list*ones(1,512);
    %bz_fluct(:,:,t_ind)=obj.bz(:,:,t_ind)-gradient(t_ind);
end
