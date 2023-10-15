clear all;
close all;

read_folder='Lx_10_Nx_256';
write_folder_list={'Lx_15_Nx_384','Lx_20_Nx_512',...
    'Lx_30_Nx_768','Lx_40_Nx_1024',...
    'Lx_80_Nx_2048','Lx_160_Nx_4096'};
Nz=64;
for write_folder_ind=1:length(write_folder_list)
    write_folder=write_folder_list{write_folder_ind};
    switch write_folder
        case 'Lx_10_Nx_256'
            Nx_write=256;
        case 'Lx_15_Nx_384'
            Nx_write=384;
        case 'Lx_20_Nx_512'
            Nx_write=512;
        case 'Lx_30_Nx_768'
            Nx_write=768;
        case 'Lx_40_Nx_1024'
            Nx_write=1024;
        case 'Lx_80_Nx_2048'
            Nx_write=2048;
        case 'Lx_160_Nx_4096'
            Nx_write=4096;
    end
    for ind=1:6
       T=h5read_complex([read_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/T'); 
       T_new=zeros(Nz,Nx_write);
       T_new(1:size(T,1),1:size(T,2))=T;
       h5write([write_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/T',T_new);

       w=h5read_complex([read_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/w'); 
       w_new=zeros(Nz,Nx_write);
       w_new(1:size(w,1),1:size(w,2))=w;
       h5write([write_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/w',w_new);

       Tz=h5read_complex([read_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/Tz'); 
       Tz_new=zeros(Nz,Nx_write);
       Tz_new(1:size(Tz,1),1:size(Tz,2))=Tz;
       h5write([write_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/Tz',Tz_new);

       wz=h5read_complex([read_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/wz'); 
       wz_new=zeros(Nz,Nx_write);
       wz_new(1:size(wz,1),1:size(wz,2))=wz;
       h5write([write_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/wz',wz_new);

       u=h5read_complex([read_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/u'); 
       u_new=zeros(Nz,Nx_write);
       u_new(1:size(u,1),1:size(u,2))=u;
       h5write([write_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/u',u_new);

       p=h5read_complex([read_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/p'); 
       p_new=zeros(Nz,Nx_write);
       p_new(1:size(p,1),1:size(p,2))=p;
       h5write([write_folder,'/X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/p',p_new);
    end
end
