function fun_10_plot_results(plot_options,raster_selected)
% ######  visual representation of the results  ######
% function to plot the results based on predefined standard plots for the
% visual representation of the most important results
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global basins_selected

breite = 20;
hoehe = 20;
schriftgroesse = 10;
visible = plot_options.plot_visibility;


% figure 1: exitcodes
if plot_options.plot_exitcodes == 1
    scale = 1;
    plot_1_spatial_exit_codes(breite*scale,hoehe*scale,schriftgroesse,visible)
end


% figure 2: spatial distribution of basins
if plot_options.plot_spatial_overview == 1
    scale = 1;
    plot_2_spatial_all_basins(breite*scale,hoehe*scale,schriftgroesse,visible,raster_selected)
end

breite = 20;
hoehe = 30;

% figure 3 & 4: single basin characteristics
if plot_options.plot_factsheet_p1 == 1 || plot_options.plot_factsheet_p2 == 1
    schriftgroesse = 10;
    scale = 1;
    for i_dam = 1:length(basins_selected)
        if plot_options.plot_factsheet_p1 == 1
            plot_3_single_basin_characteristics(i_dam,breite*scale,hoehe*scale,schriftgroesse,visible,raster_selected)
        end
        
        if plot_options.plot_factsheet_p2 == 1
            plot_4_potential_dam_characteristics_portrait(i_dam,breite*scale,hoehe*scale,schriftgroesse,visible)
        end
    end
end


% figure 5: dam characteristics
if plot_options.plot_dam_comparison == 1
    scale = 1;
    dam_ids = 1:length(basins_selected);
    plot_4_potential_dam_characteristics_portrait(dam_ids,breite*scale,hoehe*scale,schriftgroesse,visible)
end


breite = 20;
hoehe = 20;

% figure 6: basin curves
if plot_options.plot_curve_comparison == 1
    scale = 1;
    dam_ids = 1:length(basins_selected);
    plot_5_curves_depth_storage_area(dam_ids,breite*scale,hoehe*scale,schriftgroesse,visible)
end