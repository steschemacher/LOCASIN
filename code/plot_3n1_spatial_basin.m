function [sp,x_range,y_range,hContour] = plot_3n1_spatial_basin(dam_ids,schriftgroesse,fig,raster_selected,breite)
% ######  spatial plot of basin(s)  ######
% function to plot one or more dams and basins in a top view including
% background maps and elevation lines
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 23-Dez-2019

global grids grids_exclude river_points_all basins_selected


%% DATA PREPROCESSING

% extent basin and dam
if length(dam_ids)==1
    i_dam = dam_ids;
    buffer_share = .2;
    % combine coordinates
    BASIN = basins_selected(i_dam).refBASIN_depth_dam_crosssection;
    x_coords = [basins_selected(i_dam).refGK_BASIN_x(isnan(mean(BASIN,1,'omitnan'))==0),...
        basins_selected(i_dam).refGK_dam_coords(1,basins_selected(i_dam).dam_heights>0)];
    y_coords = [basins_selected(i_dam).refGK_BASIN_y(isnan(mean(BASIN,2,'omitnan'))==0),...
        basins_selected(i_dam).refGK_dam_coords(2,basins_selected(i_dam).dam_heights>0)];
    
    % get extreme values
    x_range = [min(x_coords),max(x_coords)];
    y_range = [min(y_coords),max(y_coords)];
    buffer = max(buffer_share*(x_range(2)-x_range(1)),buffer_share*(y_range(2)-y_range(1)));
    
    % define range with buffers
    x_range = [max(min(grids.x),x_range(1)-buffer),...
        min(max(grids.x),x_range(2)+buffer)];
    y_range = [max(min(grids.y),y_range(1)-buffer),...
        min(max(grids.y),y_range(2)+buffer)];
    
    % define raster data (with new extents)
    raster.x = grids.x(grids.x>=x_range(1) & grids.x<=x_range(2));
    raster.y = grids.y(grids.y>=y_range(1) & grids.y<=y_range(2));
    raster.dem = grids.dem(grids.y>=y_range(1) & grids.y<=y_range(2),...
        grids.x>=x_range(1) & grids.x<=x_range(2));
    raster.(raster_selected.name) = grids_exclude.(raster_selected.name)(grids.y>=y_range(1) & grids.y<=y_range(2),...
        grids.x>=x_range(1) & grids.x<=x_range(2));
else
    raster = grids;
    raster.(raster_selected.name) = grids_exclude.(raster_selected.name);
    x_range = [min(grids.x),max(grids.x)];
    y_range = [min(grids.y),max(grids.y)];
end

% get matrix raster data
[raster.X,raster.Y] = meshgrid(raster.x,raster.y);
raster.cellsize = grids.cellsize;

%% COLORMAPS

% schwarz-weiß (im Kreis)
anzahl = 50;
weiss_schwarz = [linspace(1,0.3,anzahl)',linspace(1,0.3,anzahl)',linspace(1,0.3,anzahl)'];
cmap_bwb = [weiss_schwarz;flipud(weiss_schwarz)];
cmap_bw = [weiss_schwarz];

% landuse
anzahl = 10;
cmap_landuse = [];
for i = 1:size(raster_selected.color,1)
   cmap_landuse = [cmap_landuse;
       repmat(raster_selected.color(i,:),anzahl,1)];
end



% caxis style
dam_height_max = ceil(max([basins_selected(dam_ids).dam_height]));
if dam_height_max <=2
    anzahl_dams = dam_height_max*2;
else
    anzahl_dams = dam_height_max;
end

caxis_range = [0,dam_height_max];
caxis_ticks = 0:dam_height_max/anzahl_dams:dam_height_max;


% dams
beige = [245,245,220]/255;
braun = [139,115,85]/255;
beige_braun = [linspace(beige(1),braun(1),anzahl_dams)',linspace(beige(2),braun(2),anzahl_dams)',linspace(beige(3),braun(3),anzahl_dams)'];
cmap_dams = [flipud(beige_braun)];


% water depths
blau1 = [105, 151, 201]/255;
blau2 = [66, 123, 184]/255;
blau3 = [51, 95, 142]/255;
blau4 = [36, 68, 101]/255;
blau14 = [linspace(blau1(1),blau4(1),anzahl_dams)',linspace(blau1(2),blau4(2),anzahl_dams)',linspace(blau1(3),blau4(3),anzahl_dams)'];
cmap_basins = [(blau14)];


% exitcodes
anzahl = 14;
gruen = [0,.8,0];
orange = [177, 108, 37]/255;
rot = [141, 36, 37]/255;
orange_rot = [linspace(orange(1),rot(1),anzahl)',linspace(orange(2),rot(2),anzahl)',linspace(orange(3),rot(3),anzahl)'];
cmap_exitcodes = [blau1];



%% background: dem
set(0,'CurrentFigure',fig)

i_subplot = 1;
sp(i_subplot) = axes(fig);%subplot(1,n_subplots,i_subplot);
hold on
[~,c] = contourf(sp(i_subplot),raster.X,raster.Y,raster.dem,'linecolor','none');
colormap(sp(i_subplot),cmap_bw)
set(sp(i_subplot),'color','none')
c.LevelStep = c.LevelStep/2;
axis equal
% define scale bar
dist_limitsX = ceil(length(raster.x)*raster.cellsize/2/100)*100;
dist_limitsY = length(raster.y)*raster.cellsize;


%% background: landuse
set(0,'CurrentFigure',fig)
i_subplot = 2;
sp(i_subplot) = axes(fig);
hold on
plot_use = imagesc(sp(i_subplot),raster.x,raster.y,raster.(raster_selected.name));
set(plot_use,'AlphaData',~isnan(raster.(raster_selected.name))*.5)
colormap(sp(i_subplot),cmap_landuse)
cb_lu = colorbar(sp(i_subplot),'Location','northoutside');
caxis([.5,length(raster_selected.legend)+.5])
set(cb_lu,'ytick',[1:length(raster_selected.legend)],...
    'yticklabel',raster_selected.legend,...
    'fontsize',schriftgroesse)
set(sp(i_subplot),'color','none')
axis equal


%% exit codes
set(0,'CurrentFigure',fig)
i_subplot = 3;
sp(i_subplot) = axes(fig);
hold on
x_coords_exitcode_all = [river_points_all.x_coord];
y_coords_exitcode_all = [river_points_all.y_coord];
x_coords_exitcode = x_coords_exitcode_all((x_coords_exitcode_all>=x_range(1) & x_coords_exitcode_all<=x_range(2)) & ...
    (y_coords_exitcode_all>=y_range(1) & y_coords_exitcode_all<=y_range(2)));
y_coords_exitcode = y_coords_exitcode_all((x_coords_exitcode_all>=x_range(1) & x_coords_exitcode_all<=x_range(2)) & ...
    (y_coords_exitcode_all>=y_range(1) & y_coords_exitcode_all<=y_range(2)));
size_x = (max(x_coords_exitcode)-min(x_coords_exitcode))/grids.cellsize+1;
size_y = (max(y_coords_exitcode)-min(y_coords_exitcode))/grids.cellsize+1;
grid_exitcode = nan(size_y,size_x);
ids_exitcode = sub2ind(size(grid_exitcode),round((y_coords_exitcode-min(y_coords_exitcode))/grids.cellsize+1),...
    round((x_coords_exitcode-min(x_coords_exitcode))/grids.cellsize+1));

grid_exitcode(ids_exitcode) = 1;
plot_exitcodes = imagesc(sp(i_subplot),min(x_coords_exitcode):grids.cellsize:max(x_coords_exitcode),...
    min(y_coords_exitcode):grids.cellsize:max(y_coords_exitcode),...
    grid_exitcode);
set(plot_exitcodes,'AlphaData',~isnan(grid_exitcode))
plot(sp(i_subplot),x_coords_exitcode,y_coords_exitcode,'.','color',cmap_exitcodes*.9,'markersize',breite*300/500)%max(length(grids.x),length(grids.y))/60)
colormap(sp(i_subplot),cmap_exitcodes*.9)
set(sp(i_subplot),'color','none')
axis equal



%% dam heights
set(0,'CurrentFigure',fig)
i_subplot = 4;
sp(i_subplot) = axes(fig);
hold on

for i_dam = dam_ids
    dam_x = min(basins_selected(i_dam).refGK_dam_coords(1,:)):grids.cellsize:max(basins_selected(i_dam).refGK_dam_coords(1,:));
    dam_y = min(basins_selected(i_dam).refGK_dam_coords(2,:)):grids.cellsize:max(basins_selected(i_dam).refGK_dam_coords(2,:));
    DAM_grid = nan(length(dam_y),length(dam_x));
    index_dam = sub2ind(size(DAM_grid),(basins_selected(i_dam).refGK_dam_coords(2,:)-min(dam_y))/grids.cellsize+1,...
        (basins_selected(i_dam).refGK_dam_coords(1,:)-min(dam_x))/grids.cellsize+1);
    DAM_grid(index_dam) = basins_selected(i_dam).dam_heights(1,:);
    plot_dam(i_dam) = imagesc(sp(i_subplot),dam_x,dam_y,DAM_grid);
    set(plot_dam(i_dam),'AlphaData',~isnan(DAM_grid))
end
colormap(sp(i_subplot),cmap_dams)
caxis(caxis_range)
cb_dam = colorbar(sp(i_subplot),'Location','westoutside','fontsize',schriftgroesse,'ytick',caxis_ticks);
ylabel(cb_dam,'dam height (m)','fontsize',schriftgroesse)
set(sp(i_subplot),'color','none')
axis equal


%% retention basins
set(0,'CurrentFigure',fig)
i_subplot = 5;
sp(i_subplot) = axes(fig);
hold on
for i_dam = dam_ids
    BASIN = basins_selected(i_dam).refBASIN_depth_dam_crosssection;
    plot_lake(i_dam) = imagesc(sp(i_subplot),basins_selected(i_dam).refGK_BASIN_x,...
        basins_selected(i_dam).refGK_BASIN_y,BASIN);
    set(plot_lake(i_dam),'AlphaData',~isnan(BASIN))  % transparent background
end
colormap(sp(i_subplot),cmap_basins)
caxis(caxis_range)

cb = colorbar(sp(i_subplot),'Location','eastoutside','fontsize',schriftgroesse,'ytick',caxis_ticks);
ylabel(cb,'water depth (m)','fontsize',schriftgroesse)
set(sp(i_subplot),'color','none')
axis equal


%% elevation lines
set(0,'CurrentFigure',fig)
i_subplot = 6;
sp(i_subplot) = axes(fig);
ax = axesm('utm');
hold on
[C,hContour] = contour(sp(i_subplot),raster.X,raster.Y,raster.dem,'ShowText','on','linecolor',[.2,.2,.2]);
clabel(C,hContour,'fontsize',schriftgroesse*.8);
hContour.LevelStep = hContour.LevelStep*2;
hContour.LabelSpacing = hContour.LabelSpacing*1.7;
set(sp(i_subplot),'color','none')
axis equal
scaleruler('units','m');
dist_colorbars = max(dist_limitsX,dist_limitsY)/35;
setm(handlem('scaleruler1'),'MajorTick', [0:dist_limitsX/5:dist_limitsX],...
    'MinorTick', 0:0,...
    'FontSize',schriftgroesse,'RulerStyle','patches','MajorTickLength',max(dist_limitsX,dist_limitsY)/80,...
    'XLoc',min(raster.x)+dist_limitsX/30,'YLoc',min(raster.y)-dist_colorbars,...
    'tickdir','down');

plot_3n2_updateContours(hContour)
addlistener(hContour, 'MarkedClean', @(h,e)plot_3n2_updateContours(hContour));

if length(dam_ids)>1
    for i_dam = dam_ids
        text(basins_selected(i_dam).x_coord , basins_selected(i_dam).y_coord , sprintf('(%d)',i_dam),'fontsize',schriftgroesse*1.1,'fontweight','bold');
    end
end
