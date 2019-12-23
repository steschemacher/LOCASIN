function plot_3_single_basin_characteristics(i_dam,breite,hoehe,schriftgroesse,visible,raster_selected)
% ######  fact sheet of one basin (page 1)  ######
% plot of the dam an basin characteristics of the best dam height
%
% functions:    fun_101_plot_spatial_one_basin
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019


global river_points_all grids basins_selected

if visible == 0
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
        'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50,'visible','off');
else
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
            'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50);    
end


% subplots 1 bis 6: spatial view of basin
position_raster = [.15,.08,.7,.4];
[sp,x_range,y_range,~] = plot_3n1_spatial_basin(i_dam,schriftgroesse,fig,raster_selected,breite*position_raster(3)*300*2.54/(fig.Position(3)));

% subplot 7: dam crosssection
set(0,'CurrentFigure',fig)
i_subplot = 7;
sp(i_subplot) = axes;
hold on
x_data = cumsum(basins_selected(i_dam).dam_axis_lengths_segments);
y_data = basins_selected(i_dam).dam_axis_heights*(-1);
if y_data(end)<0
    y_data = [y_data,0];
    x_data = [x_data,x_data(end)+x_data(end)-x_data(end-1)];
end
if y_data(1)<0
    y_data = [0,y_data];
    x_data = [0,x_data];
end
plot(x_data,y_data,'-','color',[139,90,43]/255,'linewidth',1.5);
plot(x_data([find(y_data<=0,1,'first'),find(y_data<=0,1,'last')]),[0,0],'-','color',[66, 123, 184]/255,'linewidth',1.5);
grid on
xlabel('length (m)')
ylabel('relative elevation (m)')
y_max = ceil(basins_selected(i_dam).dam_height*.5);
y_min = floor(min(y_data)*1.2);
set(gca,'ylim',[y_min,y_max]);
set(gca,'xlim',[0-.05*max(x_data),1.05*max(x_data)]);

% subplot 8: basin longitudinal section
upstream_main_river = basins_selected(i_dam).id;
distance_segment = 0;
point_dem = basins_selected(i_dam).dem;
point = basins_selected(i_dam).id;
i=0;
while point_dem < basins_selected(i_dam).dem+basins_selected(i_dam).dam_height*1.2
    i=i+1;
    points = river_points_all([river_points_all.id]==point).vorgaenger;
    if isempty(points)
        break
    end
    point = points(find([river_points_all(ismember([river_points_all.id],points)).acc]==...
        max(river_points_all(ismember([river_points_all.id],points)).acc),1,'first'));
    upstream_main_river(i) = point;
    point_dem = river_points_all(point).dem;
end
upstream_main_river_dem = [river_points_all(upstream_main_river).dem];
[river_points_all(upstream_main_river).acc];
% upstream_main_river_dem_fill = [river_points_all(upstream_main_river).dem_fill];
x_coords = [river_points_all(upstream_main_river).x_coord];
y_coords = [river_points_all(upstream_main_river).y_coord];
distance_segment = ((x_coords(2:end)-x_coords(1:end-1)).^2 + (y_coords(2:end)-y_coords(1:end-1)).^2).^(.5);

distance_cum = [0,cumsum(distance_segment)];

i_subplot = 8;
sp(i_subplot) = axes;
hold on
plot(distance_cum,upstream_main_river_dem,'-','color',[139,90,43]/255,'linewidth',1.5);
plot(([0,0]),[basins_selected(i_dam).dem,basins_selected(i_dam).dam_top_elev],'-','color',[.5,.5,.5],'linewidth',4);
plot(([0,distance_cum(find(upstream_main_river_dem<basins_selected(i_dam).dam_top_elev,1,'last'))]),...
    [basins_selected(i_dam).dam_top_elev,basins_selected(i_dam).dam_top_elev],'-','color',[66, 123, 184]/255,'linewidth',1.5);
grid on
set(gca,'xlim',[-max(distance_cum)*.05,max(distance_cum)*1.05])
set(gca,'ylim',[basins_selected(i_dam).dam_top_elev+y_min,basins_selected(i_dam).dam_top_elev+y_max])
xlabel('length (m)')
ylabel('elevation (m.a.s.l.)')

% subplot 9: water depth + area
i_subplot = 9;
sp(i_subplot) = axes;
hold on
plot(basins_selected(i_dam).ShA_dam(:,1),basins_selected(i_dam).ShA_dam(:,2),...
    '-k','linewidth',1);
plot(basins_selected(i_dam).ShA_wall(:,1),basins_selected(i_dam).ShA_wall(:,2),...
    '--k','linewidth',1);
legend('dam','wall','location','southeast','fontsize',schriftgroesse)
xlabel('basin volume (m³)        ')
ylabel('water depth (m)')
grid on

% subplot 10: water depth + volume
i_subplot = 10;
sp(i_subplot) = axes;
hold on
plot(basins_selected(i_dam).ShA_dam(:,3),basins_selected(i_dam).ShA_dam(:,2),...
    '-k','linewidth',1);
plot(basins_selected(i_dam).ShA_wall(:,3),basins_selected(i_dam).ShA_wall(:,2),...
    '--k','linewidth',1);
xlabel('basin area (m²)        ')
ylabel('water depth (m)')
grid on


%% scale and combine axes
annotation('textbox',[.05,.77,.9,.2],'String','CHARACTERISTICS AND CURVES','fontsize',schriftgroesse+2,'fontweight','bold','FitBoxToText','off')
annotation('textbox',[.05,.56,.9,.2],'String','DAM CROSS SECTION AND BASIN LONGITUDINAL SECTION','fontsize',schriftgroesse+2,'fontweight','bold','FitBoxToText','off')
annotation('textbox',[.05,.04,.9,.51],'String','SPATIAL DISTRIBUTION','fontsize',schriftgroesse+2,'fontweight','bold','FitBoxToText','off')

for i_subplot = 1:6
    set(sp(i_subplot),'XTickLabel',[],'YTickLabel',[],'fontsize',schriftgroesse)
    set(sp(i_subplot),'position',position_raster,'xlim',[x_range],'ylim',[y_range])
end
linkaxes([sp(1),sp(2),sp(3),sp(4),sp(5),sp(6)]);

set(sp(7),'position',[.12,.6,.23,.12],'fontsize',schriftgroesse)
set(sp(8),'position',[.45,.6,.45,.12],'fontsize',schriftgroesse)

set(sp(9),'position',[.45,.81,.2,.12],'fontsize',schriftgroesse)
set(sp(10),'position',[.72,.81,.2,.12],'fontsize',schriftgroesse)


dimensions = [.07,.8,.28,.15];
i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
text1 = sprintf('\n catchment area: \n dam height: \n dam axis length:\n dam volume:\n basin volume:\n specific volume:\n basin area:\n non-suited area:');

annotation('textbox',dimensions,'String',text1,'fontsize',schriftgroesse,'FitBoxToText','off','linestyle','none');
text2 = sprintf('\n %15.1f km²\n %15.1f m   \n %15.1f m   \n %15.1f m³  \n%15.1f m³  \n %15.1f mm\n %15.1f m²  \n %15.1f m²  ',...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,15),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,3),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,2),...    
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,4),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,6),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,5),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,14)*(grids.cellsize)^2);
annotation('textbox',dimensions,'String',text2,'fontsize',schriftgroesse,'FitBoxToText','off','horizontalAlignment','right','linestyle','none');

print(sprintf('single_basin_characteristics_%d.pdf',i_dam),'-dpdf','-r300')

clear fig