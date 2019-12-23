function plot_1_spatial_exit_codes(breite,hoehe,schriftgroesse,visible)
% ######  spatial distribution of the exit codes  ######
% plot to analyze the procedure of the river analysis and basin selection
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 23-Dez-2019

global grids river_points_all basins_selected

river_points_all_tmp = river_points_all;
axes_position = [.1,.1,.8,.8];
axes_xlimits = [min(grids.x),max(grids.x)];
axes_ylimits = [min(grids.y),max(grids.y)];

%%% DATA PREPROCESSING

[grids.X,grids.Y] = meshgrid(grids.x,grids.y);


%%% COLORMAPS

% schwarz-weiﬂ (im Kreis)
anzahl = 50;
weiss_schwarz = [linspace(1,0.3,anzahl)',linspace(1,0.3,anzahl)',linspace(1,0.3,anzahl)'];
cmap_bwb = [weiss_schwarz;flipud(weiss_schwarz)];
cmap_bw = [weiss_schwarz];

% landuse
anzahl = 10;
grassland = [86, 129, 45]/255;
cropland = [144, 131, 12]/255;
forest = [43, 64, 27]/255;
track = [102, 105, 99]/255;
road = [177, 108, 37]/255;
settlement = [141, 36, 37]/255;
cmap_landuse = [repmat(grassland,anzahl,1);
    repmat(cropland,anzahl,1);
    repmat(forest,anzahl,1);
    repmat(track,anzahl,1);
    repmat(road,anzahl,1);
    repmat(settlement,anzahl,1)];

% exitcodes
exit_codes_unique = unique([river_points_all_tmp.exit_code]);
anzahl_exitcodes_negative = 3;
anzahl_exitcodes_positive = 18-4;
gruen = [0,.8,0];
orange = [177, 108, 37]/255;
rot = [141, 36, 37]/255;
blau1 = [105, 151, 201]/255;
blau4 = [36, 68, 101]/255;
blau14 = [linspace(blau1(1),blau4(1),anzahl_exitcodes_negative)',linspace(blau1(2),blau4(2),anzahl_exitcodes_negative)',linspace(blau1(3),blau4(3),anzahl_exitcodes_negative)'];
blau_rot_gelb = [0,44,87,131,151,163,176,187,199,210,221,233,244,255;0,10,20,30,26,17,8,24,65,106,143,171,199,226;128,143,158,172,141,95,48,27,33,38,45,58,72,85]'/255;
cmap_exitcodes =  [0.192156866192818 0.329411774873734 0;0.258823543787003 0.496732026338577 0.0248366016894579;0.325490206480026 0.664052307605743 0.0496732033789158;0.39215686917305 0.831372559070587 0.0745098069310188;0.200000002980232 0.200000002980232 0.200000002980232;0.349999994039536 0.349999994039536 0.349999994039536;0.5 0.5 0.5;0.650000035762787 0.650000035762787 0.650000035762787;0.800000011920929 0.800000011920929 0.800000011920929;0.929411768913269 0.690196096897125 0.129411771893501;0.925490200519562 0.386274516582489 0.105882361531258;0.921568632125854 0.0823529437184334 0.0823529437184334;0.770588278770447 0.0411764718592167 0.0411764718592167;0.619607865810394 0 0;0.760784327983856 0.20392157137394 0.721568644046783;0.537254929542542 0.141176477074623 0.678431391716003;0.149019613862038 0 0.560784339904785;0 0.450980395078659 0.74117648601532;0.176470592617989 0.589542508125305 0.82745099067688;0.352941185235977 0.728104591369629 0.91372549533844;0.529411792755127 0.866666674613953 1;0.941176474094391 0.882352948188782 0.0470588244497776];


%%% FIGURE

n_subplots = 4;

if visible == 0
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
        'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50,'visible','off');
else
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
            'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50);    
end

set(0,'CurrentFigure',fig)
% DEM
i_subplot = 1;
sp(i_subplot) = axes(fig,'Position',axes_position,'xlim',axes_xlimits,'ylim',axes_ylimits);
hold on
[~,c] = contourf(sp(i_subplot),grids.X,grids.Y,grids.dem,'linecolor','none');
colormap(sp(i_subplot),cmap_bw)
c.LevelStep = c.LevelStep/2;
set(sp(i_subplot),'color','none')
axis equal
% define scale bar
dist_limitsX = ceil(length(grids.x)*grids.cellsize/2/100)*100;
dist_limitsY = length(grids.y)*grids.cellsize;


% landuse
set(0,'CurrentFigure',fig)
i_subplot = 2;
sp(i_subplot) = axes(fig,'Position',axes_position,'xlim',axes_xlimits,'ylim',axes_ylimits);
hold on
% imagesc(grids_exclude.x,grids_exclude.y,grids_exclude.landuse,'AlphaData',.5)
% % caxis([0 128])
% colormap(sp(i_subplot),cmap_landuse)
% cb_lu = colorbar('Location','northoutside');
% caxis([.5,6.5])
% set(cb_lu,'ytick',[1:6],...
%     'yticklabel',{'grassland','cropland','forest','small road','large road','settlement'},...
%     'fontsize',schriftgroesse)
set(sp(i_subplot),'color','none')
axis equal

% elevation lines
set(0,'CurrentFigure',fig)
i_subplot = 4;
sp(i_subplot) = axes(fig,'Position',axes_position,'xlim',axes_xlimits,'ylim',axes_ylimits);
axesm('utm');
hold on
[C,hContour] = contour(sp(i_subplot),grids.X,grids.Y,grids.dem,'ShowText','on','linecolor',[.2,.2,.2]);%,'AlphaData',.8)
clabel(C,hContour,'fontsize',schriftgroesse*.8);
hContour.LabelSpacing = hContour.LabelSpacing*1.7;

set(sp(i_subplot),'color','none')
axis equal

scaleruler('units','m');
setm(handlem('scaleruler1'),'MajorTick', [0:dist_limitsX/5:dist_limitsX],...
    'MinorTick', 0:0,...
    'FontSize',schriftgroesse,'RulerStyle','patches','MajorTickLength',max(dist_limitsX,dist_limitsY)/80,...
    'XLoc',min(grids.x)+dist_limitsX/30,'YLoc',min(grids.y)-max(dist_limitsX,dist_limitsY)/35,...
    'tickdir','down');

plot_3n2_updateContours(hContour)
addlistener(hContour, 'MarkedClean', @(h,e)plot_3n2_updateContours(hContour));

% exit codes

[~,index] = sortrows([river_points_all_tmp.exit_code].',-1);
river_points_all_tmp = river_points_all_tmp(index);
clear index

data_exit_code = [[river_points_all_tmp.x_coord]',[river_points_all_tmp.y_coord]',[river_points_all_tmp.exit_code]'];
data_exit_code = data_exit_code(data_exit_code(:,3)~=5,:);
i_subplot = 3;
set(0,'CurrentFigure',fig)
sp(i_subplot) = axes(fig,'Position',axes_position,'xlim',axes_xlimits,'ylim',axes_ylimits);
hold on
scatter(sp(i_subplot),data_exit_code(:,1),data_exit_code(:,2),breite*axes_position(3)*300*2.54/(fig.Position(3)),data_exit_code(:,3),'filled')
if exist('basins_selected')==1
    p = plot(sp(i_subplot),[basins_selected.x_coord],[basins_selected.y_coord],'or','markersize',10,'linewidth',2);
end
colormap(sp(i_subplot),cmap_exitcodes)
caxis([-3.5,18+.5])
set(sp(i_subplot),'color','none')
axis equal
cb = colorbar;
set(cb,'ytick',[-3:19],'fontsize',schriftgroesse)
ylabel(cb,'exit codes','fontsize',schriftgroesse)




legend(p,'selected sites','location','northeast','fontsize',schriftgroesse)

for i_subplot = 1:n_subplots
    set(sp(i_subplot),'XTickLabel',[],'YTickLabel',[])
    set(sp(i_subplot),'position',axes_position,'xlim',[min(grids.x),max(grids.x)],'ylim',[min(grids.y),max(grids.y)])
end
linkaxes([sp(1),sp(2),sp(3),sp(4)]);


print('spatial_basin_exitcodes.pdf','-dpdf','-r300')%,'-bestfit')

clear fig