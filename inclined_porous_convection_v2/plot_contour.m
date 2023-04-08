function frame=plot_contour(data,plot_config)

%%%This plot is used to plot the contour.
%%Input: data, is a cell, each cell is a structure containing x, y, z field
%%as the results we need to plot
%%plot_config, the config flag for the plotting. 


%Input:
%%data, this should be several cell. each cell should be a struct with results field 
%%Then, we will plot data{i}.x, data{i}.y, data{i}.z
%%plot_config is a struct contains the following field
%%field_all={'label_list','title_list','xlim_list','ylim_list','colormap',...
%%  'xtick_list','ytick_list','name','zlim_list','loglog','print_size','resolution','panel_num'};
%%label_list: this should be a cell, where the first one should be a
%%number, number 1 indicate we need to set the user defined label_list.
%%xlim_list: the vector of the x limit, the first is the index whether you
%%want of not.
%%ylim_list: the vecotr of the y limit, the first element is the index
%%whether you want or not. 
%%colormap: the colormap you want to specify
%%
%%zlim_list: the same as above
%%legend_list: the cell of the legend, the first element is the index
%%whether you want or not. 
%%xtick_list: the first is the index, and the following is the xtick you
%%want to set. 
%%ytick_list: the same. 
%%loglog: 1 indicate we will plot in loglog
%%panel_num: default is one. Just plot one coutour. Set as 3 to plot three
%%contour together, like convective velocity paper. 
%%print_size: the figure size that you want to print out

%%Update: 2019/08/20
%%add  'axis_2','label_list_2','xlim_list_2','ylim_list_2','xtick_list_2','ytick_list_2','loglog_2'};%%Add the default field for the plotting with double axis. Update: 2019/08/20
%%These field will be used to specify the second axis, such as showing both
%%inner unit and outer unit together...


%%resolution: the resolution to print the figure. Usually 300 is ok, some
%%journal may require 1000 resolution.

%%No output

%%---------Set default input through adding the field into the plot_config
field_all={'label_list','title_list','xlim_list','ylim_list','colormap',...
  'xtick_list','ytick_list','ztick_list','name','zlim_list','loglog','print_size','print','resolution','panel_num','marker_face_index','contour_line',...
  'axis_2','label_list_2','xlim_list_2','ylim_list_2','xtick_list_2','ytick_list_2','loglog_2','x_reverse','y_reverse','fontsize','markersize','linewidth','format','xticklabels_list','yticklabels_list','axis_equal',...
  'xtick_visible','ytick_visible','colormap_discrete','visible','streamline','colorbar','fig_num','arrow_ratio' }; %%Update 2021/02/25, add the option that I can remove the tick...
%%Add the default field for the plotting with double axis. Update: 2019/08/20
%%add the property to setup the fontsize, the default value is 40.
%%add 

field_default={{0},{0},0,0,'jet',...
    0,0,0,'test',0,[0,0],0,1,300,1,zeros(length(data)),0,...
    0,{0},0,0,0,0,[0,0],0,0,40,9,1.5,'png',{0},{0},0,...
    1,1,0,1,0,1,1,1};%%Add the default field for the plotting with double axis. Update: 2019/08/20
field_no_list=find(~isfield(plot_config,field_all)); %%fine whether there is already the field in the plot_config
for i=field_no_list
  plot_config.(field_all{i})=field_default{i};
end

close all;
if plot_config.visible
    h=figure(plot_config.fig_num);%%Plot these lines trying to reproduce figure 4 in Ahmadi et al. (2018)
else
    h=figure('Visible','Off');
end

set(gcf,'color','white');


switch plot_config.panel_num
  case {1,2,4} 
    try
        if plot_config.panel_num==1
            pcolor(data{1}.x,data{1}.y,data{1}.z); hold on;
%         elseif plot_config.panel_num==2 && 
        elseif plot_config.panel_num==4
            contour(data{1}.x,data{1}.y,data{1}.z); hold on;
        end
    catch 
       warning('data is not a matrix, plot_contour does not work'); 
    end
    if plot_config.contour_line>0
        contourf(data{1}.x,data{1}.y,data{1}.z,plot_config.contour_line); hold on;
    end
%     switch plot_config.panel_num
%       case 1 %%This plot is adding some new lines... to highlight what are these in the contour
%    
%%Update 2020/06/15, if there is more than data{2}, we also plot for the
%%quiver case.... 
     if plot_config.panel_num==2
         if plot_config.streamline==0   
            pcolor(data{1}.x,data{1}.y,data{1}.z); hold on;
            quiver(data{2}.x, data{2}.y, data{2}.u, data{2}.v,'k','Linewidth',plot_config.linewidth,'AutoScaleFactor',1.4); hold on;
            %daspect([1,1,1]); %%This is no necessary. Fix 2020/09/28, and
            %it will make some figure confusing and not consistent...
         elseif plot_config.streamline==1
             [lineobj]=streamslice(data{2}.x, data{2}.y, data{2}.u, data{2}.v,'arrow','cubic');
             %lineobj=quiver(data{2}.x, data{2}.y, data{2}.u, data{2}.v);
%              for line_ind=1:length(lineobj)
%                  lineobj(line_ind).LineWidth=3;
%              end
%              
             for line_ind=1:length(lineobj)
%                 if all(lineobj(line_ind).XData<pi)
%                     lineobj(line_ind).Color=[0,0,0];
%                     lineobj(line_ind).LineStyle='--';
%                 else
%                     lineobj(line_ind).Color=[0,0,1];
%                     lineobj(line_ind).LineStyle='-';
%                 end
                if length(lineobj(line_ind).XData)==3
                     %put the arrow as the 1/3 of the original length...
                    arrow_ratio=plot_config.arrow_ratio;
                    lineobj(line_ind).XData(1)=lineobj(line_ind).XData(2)+arrow_ratio*(lineobj(line_ind).XData(1)-lineobj(line_ind).XData(2));
                    lineobj(line_ind).XData(3)=lineobj(line_ind).XData(2)+arrow_ratio*(lineobj(line_ind).XData(3)-lineobj(line_ind).XData(2));
                    lineobj(line_ind).YData(1)=lineobj(line_ind).YData(2)+arrow_ratio*(lineobj(line_ind).YData(1)-lineobj(line_ind).YData(2));
                    lineobj(line_ind).YData(3)=lineobj(line_ind).YData(2)+arrow_ratio*(lineobj(line_ind).YData(3)-lineobj(line_ind).YData(2));

                    %delete the arrow that is in the horizontal direction
                    x1_2=lineobj(line_ind).XData(1)-lineobj(line_ind).XData(2);
                    x2_3=lineobj(line_ind).XData(2)-lineobj(line_ind).XData(3);
                    y1_2=lineobj(line_ind).YData(1)-lineobj(line_ind).YData(2);
                    y2_3=lineobj(line_ind).YData(2)-lineobj(line_ind).YData(3);
                    
                    if x1_2*x2_3<0 %|| sqrt(abs(x1_2^2+y1_2^2-x2_3^2-y2_3^2))>0.01
                        lineobj(line_ind).XData=NaN;
                        lineobj(line_ind).YData=NaN;
                    end
                    %delete the arrow that does not have the same length in
                    %the left and right 
                else
%                     lineobj(line_ind).Linewidth=3;
                end
%                 set(lineobj,'Color',[0,0,0]);
%                 set(lineobj,'MarkerSize',2);
             end
             set(lineobj,'Linewidth',plot_config.linewidth);
         elseif plot_config.streamline==2
            max_z=min(abs([max(max(data{2}.z)),min(min(data{2}.z))]));
            z_list=[-0.9*max_z, -0.7*max_z,-0.5*max_z,-0.3*max_z,-0.1*max_z,...
                0.1*max_z, 0.3*max_z, 0.5*max_z,0.7*max_z, 0.9*max_z];
            for z_ind=1:length(z_list)
                z_val=z_list(z_ind);
                line=contourc(data{2}.x,data{2}.y,data{2}.z,[z_val,z_val]); hold on;
                ind_bad=unique([find(line(1,:)>max(data{2}.x)),find(line(1,:)<min(data{2}.x)),find(line(2,:)>max(data{2}.y)), find(line(2,:)<min(data{2}.y))]);
                %ind15=unique([find(cline15(1,:)<30), find(cline15(2,:)<30),find(cline15(1,:)<5500 & cline15(2,:)>7700)]);
                line(:,ind_bad)=[];
                
                clear line_list;
                line_list_ind=1;
                line_list_sub_ind=1;
                for line_ind=1:length(line)-1
                    line_list{line_list_ind}(:,line_list_sub_ind)=line(:,line_ind);
                    if norm(line(:,line_ind+1)-line(:,line_ind))>0.1
                        line_list_ind=line_list_ind+1;
                        line_list_sub_ind=1;
                    else
                       line_list_sub_ind=line_list_sub_ind+1;
                    end
                end
                line_ind=line_ind+1;
                line_list{line_list_ind}(:,line_list_sub_ind)=line(:,line_ind);
                
                if z_val==0
                   %plot(line(1,:),line(2,:),'r-.'); 
                elseif z_val>0
                    for line_list_ind=1:length(line_list)
                       plot(line_list{line_list_ind}(1,:),line_list{line_list_ind}(2,:),'k-','Linewidth',plot_config.linewidth); 
                    end
%                     ind1=find(line(1,:)<pi);
%                     line1=line(:,ind1);
%                     ind2=find(line(1,:)>pi);
%                     line2=line(:,ind2);
%                     plot(line1(1,:),line1(2,:),'k-','Linewidth',plot_config.linewidth);
%                     plot(line2(1,:),line2(2,:),'k-','Linewidth',plot_config.linewidth);

                elseif z_val<0
                    for line_list_ind=1:length(line_list)
                       plot(line_list{line_list_ind}(1,:),line_list{line_list_ind}(2,:),'b--','Linewidth',plot_config.linewidth); 
                    end
%                     ind1=find(line(1,:)<pi);
%                     line1=line(:,ind1);
%                     ind2=find(line(1,:)>pi);
%                     line2=line(:,ind2);
%                     plot(line1(1,:),line1(2,:),'b--','Linewidth',plot_config.linewidth);
%                     plot(line2(1,:),line2(2,:),'b--','Linewidth',plot_config.linewidth);

                end
            end
            grid on; %set(gca,'gridlinestyle','Linewidth',1.5);
            grid minor;  %%set up the minor grid.
            ax=gca;
            ax.GridAlpha=0.55;
            ax.MinorGridAlpha=0.35;
         end
         
         
     end
    if length(data)>plot_config.panel_num
       for lin_ind=(plot_config.panel_num+1):length(data)
        %%Note this marker size is modified... Which is not consistent with
        %%the plot_line
        tmp=plot(data{lin_ind}.x,data{lin_ind}.y,...
            plot_config.user_color_style_marker_list{lin_ind-plot_config.panel_num},'Linewidth',plot_config.linewidth,'Markersize',plot_config.markersize); hold on;
        if plot_config.marker_face_index(lin_ind-plot_config.panel_num)
            set(tmp,'MarkerFaceColor',plot_config.marker_face_color{lin_ind-plot_config.panel_num});
            set(tmp,'MarkerEdgeColor',plot_config.marker_face_color{lin_ind-plot_config.panel_num});
        end
       end
    end
        %%These are fixed 2019/07/10. It will show the tick on top of the contour
        %%plots... If we comment the grid off, the grid will also show on top of
        %%the contour.
        set(gca, 'Layer', 'top');
        grid off;
     
    
    if plot_config.label_list{1} %%Set up the label
      try
        xlabel(plot_config.label_list{2},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
        ylabel(plot_config.label_list{3},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
      catch 
        disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
        xlabel(plot_config.label_list{2},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
        ylabel(plot_config.label_list{3},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
      end
    end
    if plot_config.loglog(1)
      set(gca,'xscale','log');
    end
    if plot_config.loglog(2)
      set(gca,'yscale','log');
    end 
    if plot_config.xlim_list(1)
      set(gca,'xlim',[plot_config.xlim_list(2),plot_config.xlim_list(3)]);
    end
    if plot_config.ylim_list(1)
      set(gca,'ylim',[plot_config.ylim_list(2),plot_config.ylim_list(3)]);
    end
    if plot_config.xtick_list(1) %%Set up x, y tick
      set(gca,'xtick',plot_config.xtick_list(2:end));
    end
    if plot_config.ytick_list(1)
      set(gca,'ytick',plot_config.ytick_list(2:end));
    end
    if plot_config.xticklabels_list{1}
        set(gca,'xticklabels',plot_config.xticklabels_list(2:end))
    end
    if plot_config.yticklabels_list{1}
         set(gca,'yticklabels',plot_config.yticklabels_list(2:end))
    end
%     if plot_config.zticklabels_list{1}
%          set(gca,'cticklabels',plot_config.zticklabels_list(2:end))
%     end
    if plot_config.title_list{1} %%Set up title
       try 
        title(plot_config.title_list{2},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize)
       catch 
        disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
        title(plot_config.title_list{2},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize)
       end
    end
    if plot_config.zlim_list(1) %%set up x, y, z limit
      caxis([plot_config.zlim_list(2),plot_config.zlim_list(3)]);
    end
    try
       set(gca, 'TickLabelInterpreter','Latex', 'Fontsize', plot_config.fontsize);
    catch
       disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
       set(gca, 'Fontsize', plot_config.fontsize);
    end
    shading interp;
    if plot_config.colorbar
        try 
            cb=colorbar('TickLabelInterpreter','Latex', 'Fontsize', plot_config.fontsize);
            set(gca, 'TickLabelInterpreter','Latex', 'Fontsize', plot_config.fontsize);
        catch
            disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
            colorbar( 'Fontsize', plot_config.fontsize);
            set(gca, 'Fontsize', plot_config.fontsize);
        end

        if plot_config.ztick_list(1)
            cb.Ticks=plot_config.ztick_list(2:end);
        end
    
    end
    
    cb.Ruler.Exponent = 0;

%     set(h, 'Position', [.05 .15 .9 .05]);

    
    
    if plot_config.print_size(1)
      set(gcf,'outerposition',[1,1,plot_config.print_size(2:3)]);
    end
%     set(cb, 'YTickLabel', cellstr(num2str(reshape(get(cb, 'YTick'),[],1),'%0.3f')) )
%      set(cb, 'YTickLabel', cellstr(num2str(reshape(get(cb, 'YTick'),[],1),'%6s')));

    %%These seems to be very difficult to work in a nice way....
    if plot_config.axis_2 %%Set up the double axis... y label on the right, x label on the top
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');
        if plot_config.label_list_2{1} %%Set up the label
          try
            xlabel(plot_config.label_list_2{2},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
            ylabel(plot_config.label_list_2{3},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
          catch 
            disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
            xlabel(plot_config.label_list_2{2},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
            ylabel(plot_config.label_list_2{3},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
          end
        end
        if plot_config.loglog_2(1)
          set(ax2,'xscale','log');
        end
        if plot_config.loglog_2(2)
          set(ax2,'yscale','log');
        end 
        if plot_config.xlim_list_2(1)
          set(ax2,'xlim',[plot_config.xlim_list_2(2),plot_config.xlim_list_2(3)]);
        end
        if plot_config.ylim_list_2(1)
          set(ax2,'ylim',[plot_config.ylim_list_2(2),plot_config.ylim_list_2(3)]);
        end
        if plot_config.xtick_list_2(1) %%Set up x, y tick
          set(ax2,'xtick',plot_config.xtick_list_2(2:end));
        end
        if plot_config.ytick_list_2(1)
          set(ax2,'ytick',plot_config.ytick_list_2(2:end));
        end
        try
           set(gca, 'TickLabelInterpreter','Latex', 'Fontsize', plot_config.fontsize);
        catch
           disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
           set(gca, 'Fontsize', plot_config.fontsize);
        end
        try 
            %cb=colorbar('TickLabelInterpreter','Latex', 'Fontsize', plot_config.fontsize);
            set(ax2, 'TickLabelInterpreter','Latex', 'Fontsize', plot_config.fontsize);
        catch
            disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
            %colorbar( 'Fontsize', plot_config.fontsize);
            set(gca, 'Fontsize', plot_config.fontsize);
        end
%         pos = get(cb,'Position');
%         set(cb,'Location','south');
%         11 
%         set(gcf,'Position');
    end
    
%     if plot_config.print_size(1)
%       set(gcf,'outerposition',[1,1,plot_config.print_size(2:3)]);
%     end
  case 3 %%Plot three panel together
    set(gcf,'outerposition',[1,1,1650,830]);
    plot_config.print_size=0;
    set(gcf,'color','white');
    clf();
    hold on;
    % Get the width and height of the figure
    lbwh = get(1, 'position');
    figw = lbwh(3);
    figh = lbwh(4);
    % Number of rows and columns of axes
    ncols = 3;
    nrows = 1;
    % w and h of each axis in normalized units
    axisw = (1 / ncols) * 0.85;
%     axisw = (1 / ncols) * 0.75;
    axish = (1 / nrows) * 0.555;
%     axish = (1 / nrows) * 0.455;

    indexplot=1;  
    
    for i=1:plot_config.panel_num
      row = floor( indexplot/(ncols+1))+1;
      col = mod( indexplot-1, ncols )+1;
      axisl = (axisw+0.01) * (col-1)+0.12;
%       axisl = (axisw+0.01) * (col-1)+0.22;
      axisb = (axish+0.01) * (row-1)+0.35;
      subplot('Position',[axisl,axisb,axisw,axish]);
      pcolor(data{i}.x,data{i}.y,data{i}.z); hold on;
      %%These are fixed 2019/07/10. It will show the tick on top of the contour
        %%plots... If we comment the grid off, the grid will also show on top of
        %%the contour.
 
      %%Plot the line over the contour if any 
      if length(data)>3
        for lin_ind=4:length(data)
        plot(data{lin_ind}.x,data{lin_ind}.y,...
            plot_config.user_color_style_marker_list{lin_ind-plot_config.panel_num},'Linewidth',plot_config.linewidth,'Markersize',plot_config.markersize); hold on; %original size is 12
        end
      end
      
      if (indexplot==1) %%set up the y tick and y label only for the first one
          if plot_config.ytick_list(1)
            set(gca,'ytick',plot_config.ytick_list(2:end));
          end
          if plot_config.label_list{1}
            try
              ylabel(plot_config.label_list{3},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
            catch
              disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
              ylabel(plot_config.label_list{3},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
            end
            
          end
      else
          set(gca,'ytick',[]);
      end
      %%Set up x tick and x label
      if plot_config.xtick_list(1)
        set(gca,'xtick',plot_config.xtick_list(2:end));
      end
      if plot_config.label_list{1} %Set up x label
        try
          xlabel(plot_config.label_list{2},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
        catch
          disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
          xlabel(plot_config.label_list{2},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize);
        end
      end
      
      %%Set up x, y, z lim
      if plot_config.xlim_list(1)
        xlim([plot_config.xlim_list(2),plot_config.xlim_list(3)]);
      end
      if plot_config.ylim_list(1)
          ylim([plot_config.ylim_list(2),plot_config.ylim_list(3)]);
      end
      if plot_config.zlim_list(1)
          caxis([plot_config.zlim_list(2),plot_config.zlim_list(3)]);
      end
      
      %%Set up title of subfigures
      if plot_config.title_list{1}
          try
            title(plot_config.title_list{indexplot+1},'Interpreter','Latex','FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize)
          catch
            disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
            title(plot_config.title_list{indexplot+1},'FontWeight','bold','Fontname', 'Times New Roman', 'Fontsize', plot_config.fontsize)
          end
      end
      if plot_config.loglog
        set(gca,'xscale','log');
        set(gca,'yscale','log');
      end
      try
        set(gca, 'TickLabelInterpreter','Latex', 'Fontsize', plot_config.fontsize);
      catch
        disp('Latex interpreter is not supported on this machine. Use Time New Roman font instead.');
        set(gca, 'Fontsize', plot_config.fontsize);
      end
      
        
      indexplot=indexplot+1;
      shading interp;
      
      set(gca, 'Layer', 'top');
       grid off;
      
    end
    cb=colorbar('south','Position',[0.118,0.1,0.858,0.05],'TickLabelInterpreter','Latex', 'Fontsize', plot_config.fontsize);
    if plot_config.ztick_list(1)
        cb.Ticks=plot_config.ztick_list(2:end);
    end
    %set(cb, 'YTickLabel', cellstr(num2str(reshape(get(cb, 'YTick'),[],1),'%0.3f')) )

    %%Update 2020/09/29, take this out of switch panel number
    if plot_config.print_size(1)
      set(gcf,'outerposition',[1,1,plot_config.print_size(2:3)]);
    end
    otherwise 
    error('Wrong panel_num');
        
end


% 
% grid on; %set(gca,'gridlinestyle','Linewidth',1.5);
% ax=gca;
% ax.GridAlpha=0.6;
% grid off;

if plot_config.x_reverse
    set ( gca, 'xdir', 'reverse');
end
if plot_config.y_reverse
    set ( gca, 'ydir', 'reverse');
end

%%Update 2020/12/17, add the setup to make axis in the equal spacing
if plot_config.axis_equal
    axis equal;
end

if plot_config.colormap_discrete==0
    switch plot_config.colormap %%Set up colorbar
      case 'jet'
          colormap(jet);
      case 'white_zero'
    %     load('./post_processing/MyColormaps.mat');
         load('MyColormaps.mat'); 
         %%Update 2019/9/10, fix the bug when then main function calling this
         %%is not in the main folder... The above may fail.

        colormap(mycmap);
      case 'redblue'
        addpath('./post_processing');
        colormap(redblue);
      case 'bluewhitered' 
        addpath('./post_processing');%%Note 2020/03/18, this colormap will enforce zero as white color, no matter 0 is the middle point or not
        colormap(bluewhitered); shading interp;
      case 'gray'
        shading interp;  colormap(gray);
        case 'hot' %%add the option for the hot color
        shading interp; colormap(hot);
        case 'default'

      otherwise 
        error('Wrong color map');
    end
    
else%%This is setting the discrete colormap numbers.
      switch plot_config.colormap %%Set up colorbar
      case 'jet'
         colormap(jet(plot_config.colormap_discrete));
      case 'white_zero'
    %     load('./post_processing/MyColormaps.mat');
         load('MyColormaps.mat'); 
         %%Update 2019/9/10, fix the bug when then main function calling this
         %%is not in the main folder... The above may fail.

        colormap(mycmap(plot_config.colormap_discrete));
      case 'redblue'
        addpath('./post_processing');
        colormap(redblue(plot_config.colormap_discrete));
      case 'bluewhitered' 
        addpath('./post_processing');%%Note 2020/03/18, this colormap will enforce zero as white color, no matter 0 is the middle point or not
        colormap(bluewhitered(plot_config.colormap_discrete)); shading interp;
      case 'gray'
        shading interp;  colormap(gray(plot_config.colormap_discrete));
      case 'hot' %%add the option for the hot color
        shading interp; colormap(hot(plot_config.colormap_discrete));
      case 'default'
        shading interp;
        colormap(parula(plot_config.colormap_discrete));
      otherwise 
        error('Wrong color map');
    end
end

if plot_config.print_size(1)
      set(gcf,'outerposition',[1,1,plot_config.print_size(2:3)]);
end

if ~plot_config.xtick_visible
    set(gca,'Xticklabel',[]);
end

if ~plot_config.ytick_visible
    set(gca,'Yticklabel',[]);
end



if plot_config.print
    set(gcf,'PaperPositionMode', 'auto'); 
    switch plot_config.format
        case 'png'
             print(gcf,'-dpng',strcat('-r',num2str(plot_config.resolution)),plot_config.name);
        case 'eps'
             print(gcf,'-depsc',strcat('-r',num2str(plot_config.resolution)),plot_config.name);
        case 'tiff'
             print(gcf,'-depsc','-tiff',strcat('-r',num2str(plot_config.resolution)),plot_config.name);
    end
end

frame=getframe(h);

end

