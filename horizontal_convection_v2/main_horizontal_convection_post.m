clear all;
close all;
clc;

%slurm_num={'20230408221941'};
slurm_num={'14260113'}; %change this number as the output in your file.
flag.print=0;
flag.visible=0;
flag.video=1;
flag.no_ylabel=0;
for slurm_ind=1:length(slurm_num)
    if length(slurm_num{slurm_ind})==14
        h5_name=['D:\Data\dedalus_example\inclined_porous_convection_v2\analysis_',...
            slurm_num{slurm_ind},'\analysis_s1.h5'];
    elseif length(slurm_num{slurm_ind})==8
        h5_name=['D:\Data\dedalus_example\inclined_porous_convection_v2\dedalus_',...
            slurm_num{slurm_ind},'\analysis\analysis_s1.h5'];
    end
    dedalus_post_my{slurm_ind}=horizontal_convection_post(h5_name,flag);
%     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('u',[],[],[100]);
%     dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('w',[],[],[100]);
    dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.snapshot('T_tot',[],[],[]);
    %dedalus_post_my{slurm_ind}=dedalus_post_my{slurm_ind}.spectrum_average('T');

end
