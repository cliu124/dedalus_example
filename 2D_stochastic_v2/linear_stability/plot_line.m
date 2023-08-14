function frame=plot_line(data,plot_config)
%plot_line Summary of this function goes here
%   Detailed explanation goes here

%%Chang Liu
%%2018/12/15
%%Updated: 2019/03/04, I add the input to specify the printed figure size.
%%Update: 2019/03/15, simplify the input as the struct of data and
%%plot_config
%Input:
%%data, this should be several cell. each cell contain a two row matrix. 
%%Then, wmine will plot data{i}(1,:), data{i}(2,:).
%%plot_config is a struct contains the following field
%% label_list,title_list,xlim_list,ylim_list,legend_list,xtick_list,ytick_list,name,Markerindex,user_color_style_marker_list,print_size,resolution
%%label: this should be a one by two vector, it give the string of the x
%%label, ylabel and
%%xlim_list: the list of the x limit, the first is the index whether you
%%want of not.
%%ylim_list: the list of the y limit, the first element is the index
%%whether you want or not. 
%%legend_list: the list of the legend, the first element is the index
%%whether you want or not. 
%%xtick_list: the first is the index, and the following is the xtick you
%%want to set. 
%%ytick_list: the same. 


%%Markerindex: index 1, only use the line without marker
%%index 2, use both the line and the marker
%%index 3, use the last input, user defined input. 

%%user_color_style_marker_list: if not, just input 0. 

%%print_size: the figure size that you want to print out
%%

%%resolution

%%No output, just print the figure at current file folder.

%%Sample usage
% % 
% % for i=1:length(a)
% % plot_config.Markerindex=3;   
% % plot_config.user_color_style_marker_list={'k-','b--','rx'}; plot_config.resolution=300;
% % data{1}=[time;state_uncon(i,:)]; data{2}=[time;state_con(i,:)];
% % plot_config.label_list={1,'$t$',strcat('$a_',num2str(i),'(t)$')};
% % plot_config.title_list={0};
% % plot_config.legend_list={0,'No feedback','With feedback'};
% % plot_config.xlim_list=[1,0,1000]; plot_config.ylim_list=[1,-1.5,1.5]; plot_config.xtick_list=[1,0,200,400,600,800,1000]; plot_config.ytick_list=[1,-1.5,-1,-0.5,0,0.5,1,1.5];
% % plot_config.name=strcat('passive_NS_9D_Galerkin_mean1_9_state',num2str(i),'_20190210.png');
% % plot_config.print_size=[1,1200,1100];
% % plot_line(data,plot_config);
% % end


%%---------Set default input:
field_all={'label_list','title_list','xlim_list','ylim_list','legend_list',...
  'xtick_list','ytick_list','name','Markerindex','user_color_style_marker_list'...
  ,'print_size','print','resolution','loglog','marker_face_index','fontsize','fontsize_legend','markersize','linewidth','xticklabels_list','yticklabels_list','visible','legend_index','RGB'}; %%add the linewidth flag and default value is 1.5
field_default={{0},{0},0,0,{0},0,0,'test',1,0,0,1,300,[0,0],zeros(length(data)),40,40,12,1.5,{0},{0},1,0,{0}};
field_no_list=find(~isfield(plot_config,field_all)); %%fine whether there is already the field in the plot_config
for i=field_no_list
  plot_config.(field_all{i})=field_default{i};
end

%%-------------
%%These are default color style and marker style
color_style_list={'k-','b--','r-.','m:'};
color_style_marker_list={'k-o','b--+','r-.*','m:x'};

switch plot_config.Markerindex
  case 1
    this_color=color_style_list;
  case 2
    this_color=color_style_marker_list;
  case 3
    this_color=plot_config.user_color_style_marker_list;
  otherwise
    disp('Please input correct marker index');
end
close all;

if plot_config.visible
    h=figure;%%Plot these lines trying to reproduce figure 4 in Ahmadi et al. (2018)
else
    h=figure('Visible','Off');
end

%h=figure(1)%%Plot these lines trying to reproduce figure 4 in Ahmadi et al. (2018)
set(gcf,'color','white');

%%The data are stored as the column vector
for i=1:length(data)
  tmp=plot(data{i}.x,data{i}.y,this_color{i},'Linewidth',plot_config.linewidth,'Markersize',plot_config.markersize); hold on;
  if plot_config.RGB{1}%%Update 2022/03/01, add the RGB setting... 
      tmp.Color=plot_config.RGB{i+1};
  end
  plot_handle.(['index',num2str(i)])=tmp;
  %%Marker face color, to plot the filled marker
  if plot_config.marker_face_index(i)
        set(tmp,'MarkerFaceColor',plot_config.marker_face_color{i});
        set(tmp,'MarkerEdgeColor',plot_config.marker_face_color{i});
  end
end
grid on; %set(gca,'gridlinestyle','Linewidth',1.5);
grid minor;  %%set up the minor grid.
ax=gca;
ax.GridAlpha=0.55;
ax.MinorGridAlpha=0.35;

if plot_config.xlim_list(1)
  xlim([plot_config.xlim_list(2),plot_config.xlim_list(3)]);
end
if plot_config.ylim_list(1)
  ylim([plot_config.ylim_list(2),plot_config.ylim_list(3)]);
end
if plot_config.xtick_list(1)
  set(gca,'xtick',plot_config.xtick_list(2:end));
end
if plot_config.ytick_list(1)
  set(gca,'ytick',plot_config.ytick_list(2:end));
end
if plot_config.loglog(1)
      set(gca,'xscale','log');
end
if plot_config.loglog(2)
      set(gca,'yscale','log');
end 
if plot_config.label_list{1}
  try
    xlabel(plot_config.label_list{2},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
    ylabel(plot_config.label_list{3},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
  catch
    disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
    xlabel(plot_config.label_list{2},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
    ylabel(plot_config.label_list{3},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
  end
end
if plot_config.title_list{1}
    title(plot_config.title_list{2},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
end

if plot_config.legend_list{1}
  if plot_config.legend_index(1)==0
     plot_config.legend_index=1:(length(plot_config.legend_list)-1);
  end
  for ind=1:length(plot_config.legend_index) 
      legend_index_tmp(ind)=plot_handle.(['index',num2str(plot_config.legend_index(ind))]);
  end
  lg=legend(legend_index_tmp,plot_config.legend_list(2:end));
  try
  set(lg,'Interpreter','Latex','Location','best','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize_legend);
  catch 
    disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
    set(lg,'Location','best','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize_legend);
  end
end
if plot_config.print_size(1)
  set(gcf,'outerposition',[1,1,plot_config.print_size(2:3)]);
else
  set(gcf,'outerposition',[1,1,1600,1200]);
end

try
    set(gca, 'TickLabelInterpreter','Latex', 'Fontsize', plot_config.fontsize);
catch
    disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
    set(gca,'Fontsize', plot_config.fontsize);
end
if plot_config.print
    set(gcf,'PaperPositionMode', 'auto');
    print(gcf,'-dpng',strcat('-r',num2str(plot_config.resolution)),plot_config.name);
end

frame=getframe(h);

end

