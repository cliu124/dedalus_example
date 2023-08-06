clear all;
close all;

for ind=1:5
   T=h5read_complex(['X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/T'); 
   h5write(['X-',num2str(ind),'_checkpoint_s1.h5'],'/tasks/T',-fliplr(flipud(T)));
   
   w=h5read_complex(['X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/w'); 
   h5write(['X-',num2str(ind),'_checkpoint_s1.h5'],'/tasks/w',-fliplr(flipud(w)));
   
   Tz=h5read_complex(['X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/Tz'); 
   h5write(['X-',num2str(ind),'_checkpoint_s1.h5'],'/tasks/Tz',fliplr(flipud(Tz)));
   
   wz=h5read_complex(['X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/wz'); 
   h5write(['X-',num2str(ind),'_checkpoint_s1.h5'],'/tasks/wz',fliplr(flipud(wz)));
   
   u=h5read_complex(['X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/u'); 
   h5write(['X-',num2str(ind),'_checkpoint_s1.h5'],'/tasks/u',fliplr(flipud(u)));
   
   p=h5read_complex(['X',num2str(ind),'_checkpoint_s1.h5'],'/tasks/p'); 
   h5write(['X-',num2str(ind),'_checkpoint_s1.h5'],'/tasks/p',fliplr(flipud(p)));
   
end