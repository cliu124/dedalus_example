clear all;
close all;
clc;

slurm_num={'20230406131202'};
flag.print=1;
flag.visible=1;
flag.video=1;
flag.no_ylabel=0;
for slurm_ind=1:length(slurm_num)
    h5_name=['D:\Data\dedalus_example\2D_stochastic_v2\snapshots_',...
        slurm_num{slurm_ind},...
        '\snapshots_s1.h5'];
    stochastic_post_my{slurm_ind}=stochastic_post(h5_name,flag);
    %var=h5read(h5_name,'/tasks/forcing_var_x_coeff');
    %var=h5read(h5_name,'/tasks/forcing_var_x');
    
end
