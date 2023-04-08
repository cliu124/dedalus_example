function plot_video(frame_list,plot_config)

%%Update 2021/09/01
%%plot video..
video=VideoWriter(plot_config.name,'MPEG-4');
video.FrameRate=10;
open(video);
writeVideo(video,frame_list);
close(video);
end

