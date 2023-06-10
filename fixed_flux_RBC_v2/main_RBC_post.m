clear all;
close all;
clc;

slurm_num={'20230406211615'};
flag.print=0;
flag.visible=0;
flag.video=0;
flag.no_ylabel=0;
for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus_example\fixed_flux_RBC_v2\analysis_',...
        slurm_num{slurm_ind},'\analysis_s1.h5'];
    dedalus_post_my{slurm_ind}=RBC_post(h5_name,flag);
    %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('u',[],[],[10:10:100]);
    %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('w',[],[],[10:10:100]);
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('b',[],[],[10:10:100]);
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('b');

end
