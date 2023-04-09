clear all;
close all;
clc;

slurm_num={'20230408221941'};
slurm_num={'20230408224303'};
flag.print=1;
flag.visible=1;
flag.video=0;
flag.no_ylabel=0;
for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus_example\inclined_porous_convection_v2\analysis_',...
        slurm_num{slurm_ind},'\analysis_s1.h5'];
    dedalus_post_my{slurm_ind}=post_inclined_porous_convection(h5_name,flag);
%     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('u',[],[],[100]);
%     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('w',[],[],[100]);
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T_tot',[],[],[100]);
    %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('T');

end
