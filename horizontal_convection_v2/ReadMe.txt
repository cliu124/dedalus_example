This assumes you are running on the Ubuntu virtual machine of windows system, and dedalus2 has been installed following 
https://dedalus-project.readthedocs.io/en/v2_master/pages/installation.html

1. On the ubuntu terminal, type

sh run.sh

will then run the main code. This will also automatically copy the data to windows system that is D:/Data/dedalus_example/inclined_porous_convection_v2/

2. If you want to stop the run, press ctrl+c several times. Then, you can continue post-processing by typing in the terminal

sh post.sh 

3. For the post-processing to generate the video, run the MATLAB code (this can be down in windows system)

main_inclined_porous_convection.m

This code will generate the picture and video within the folder of analysis_yyyyMMddhhmmss

4. For the control parameter of this problem, you can modify that in the line 48-69 of the code

inclined_porous_convection.py

5. To modify the initial condition that continue from, modify the second line in run.sh

by changing 'X1_checkpoint_s1.h5' into other h5 file names within this folder. 

'X1_checkpoint_s1.h5' localized structure with one pulse
'X2_checkpoint_s1.h5' localized structure with two pulses
'X3_checkpoint_s1.h5' localized structure with three pulses
'X4_checkpoint_s1.h5' localized structure with four pulses
'X5_checkpoint_s1.h5' steady convection roll with five rolls