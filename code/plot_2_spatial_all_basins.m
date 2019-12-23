function plot_2_spatial_all_basins(breite,hoehe,schriftgroesse,visible,raster_selected)
% ######  spatial distribution of dams and basins  ######
% function to plot locations of the dams and basins of the selected basin
% combination including dam heights and water depths
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global basins_selected 
position_raster = [.1,.1,.8,.8];

if visible == 0
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
        'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50,'visible','off');
else
    fig = figure('color',[1,1,1],'rend','painters','paperunits','centimeters','papersize',[breite,hoehe],...
            'paperpositionmode','manual','paperPosition',[0 0 breite hoehe],'position',[1 1 breite hoehe]*50);    
end

set(0,'CurrentFigure',fig)
n_subplots = 6;
dam_ids = 1:length(basins_selected);
% subplots 1 bis 6: spatial view of basin
[sp,x_range,y_range,~] = plot_3n1_spatial_basin(dam_ids,schriftgroesse,fig,raster_selected,breite*position_raster(3)*300*2.54/(fig.Position(3)));

for i_subplot = 1:6
    set(sp(i_subplot),'XTickLabel',[],'YTickLabel',[])
    set(sp(i_subplot),'position',position_raster,'xlim',[x_range],'ylim',[y_range])
end

linkaxes([sp(1),sp(2),sp(3),sp(4),sp(5),sp(6)]);


print('spatial_all_basins.pdf','-dpdf','-r300')%,'-bestfit') 

clear fig sp