clear all;
close all;
clc;

slurm_num={'20230406165535'};
flag.print=1;
flag.visible=1;
flag.video=0;
flag.no_ylabel=0;
for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus_example\double_diffusive_settling_v2\analysis_',...
        slurm_num{slurm_ind},...
        '\analysis_s1.h5'];
    dedalus_post_my{slurm_ind}=settling_post(h5_name,flag);

    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('u',[],[],[10:10:100]);
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('w',[],[],[10:10:100]);
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T',[],[],[10:10:100]);
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('C',[],[],[10:10:100]);

    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('C');
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('T');

    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_TKE_average();


end
