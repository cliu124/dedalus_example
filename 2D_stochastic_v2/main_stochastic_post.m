clear all;
close all;
clc;

slurm_num={'20230406135106'};
flag.print=1;
flag.visible=1;
flag.video=0;
flag.no_ylabel=0;
for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus_example\2D_stochastic_v2\analysis_',...
        slurm_num{slurm_ind},...
        '\analysis_s1.h5'];
    stochastic_post_my{slurm_ind}=stochastic_post(h5_name,flag);
    stochastic_post_my{slurm_ind}=stochastic_post_my{slurm_ind}.snapshot('u',[],[],[10:10:100]);
    stochastic_post_my{slurm_ind}=stochastic_post_my{slurm_ind}.snapshot('v',[],[],[10:10:100]);
    stochastic_post_my{slurm_ind}=stochastic_post_my{slurm_ind}.snapshot('forcing_var_x',[],[],[10:10:100]);
    stochastic_post_my{slurm_ind}=stochastic_post_my{slurm_ind}.snapshot('forcing_var_y',[],[],[10:10:100]);

    stochastic_post_my{slurm_ind}=stochastic_post_my{slurm_ind}.spectrum_snapshot('u');
    stochastic_post_my{slurm_ind}=stochastic_post_my{slurm_ind}.spectrum_average('u');
    stochastic_post_my{slurm_ind}=stochastic_post_my{slurm_ind}.spectrum_snapshot('forcing_var_x');

    stochastic_post_my{slurm_ind}=stochastic_post_my{slurm_ind}.spectrum_average('forcing_var_x');

    stochastic_post_my{slurm_ind}=stochastic_post_my{slurm_ind}.spectrum_TKE_average();


    %var=h5read(h5_name,'/tasks/forcing_var_x_coeff');
    %var=h5read(h5_name,'/tasks/forcing_var_x');
    
end