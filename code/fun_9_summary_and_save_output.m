function fun_9_summary_and_save_output(save_options)
% ######  save results  ######
% function to summary the output and to save the results based on the user
% definitions (selected detail of the output)
%
% functions:    fun_91_save_basins_as_ascii.m
%               fun_92_save_basins_characteristics_as_xlsx
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global path_result refGK_LLcornerGLOBAL dam_height_max grids river_points_all
global dam_points basins_selected
global  grids_exclude info_exclude_dam info_exclude_basin 

cd(path_result)

%%% used grids
% grids, which were imported and used for the analysis
if save_options.save_grids == 1    
    save('input_grids_used.mat','grids','grids_exclude',...
        'info_exclude_dam','info_exclude_basin','-v7.3')
end


%%% river points
%   basic information for all river points
if save_options.save_river_points == 1
    fields = setdiff(fieldnames(river_points_all),{'id','x_coord','y_coord','exit_code','dem','acc','dem_fill','vorgaenger','nachfolger'});
    river_points_all = rmfield(river_points_all,fields);
    save('river_points.mat','river_points_all','-v7.3')
end


%%% dam_points
% parameters for all possible dam sites
if save_options.save_dam_points == 1
    fields = {'dam_top_elev','dam_height','dam_heights','dam_elevations','dam_volumes',...
        'dam_axis_lengths','dam_axis_heights','basin_volumes','basin_areas',...
        'crit1_damvol_per_storage', 'crit2_area_per_storage', ...
        'crit3_share_wellsuited_cells', 'crit4_share_notsuited_cells'};
    for i_basins = 1:length(dam_points)
        for i_field = 1:length(fields)
            dam_points(i_basins).(fields{i_field}) = dam_points(i_basins).(fields{i_field})(dam_points(i_basins).dam_row_selected,:);
        end
    end
    fields = setdiff(fieldnames(dam_points),...
        {'id', 'x_coord', 'y_coord', 'refGK_dam_coords', ...
        'dem','acc','dam_top_elev', 'dam_height','dam_axis_lengths',...
        'dam_axis_lengths_segments','dam_axis_heights',...
        'refGK_BASIN_x','refGK_BASIN_y',...
        'refBASIN_depth_dam_crosssection', 'refBASIN_depth_wall_crosssection',...
        'dam_heights','dam_elevations','dam_volumes','basin_volumes','basin_areas',...
        'crit1_damvol_per_storage', 'crit2_area_per_storage', ...
        'crit3_share_wellsuited_cells', 'crit4_share_notsuited_cells','objective_function',...
        'ShA_dam', 'ShA_wall','curve_potential_dam_basin'});
    dam_points = rmfield(dam_points,fields);
    save('dam_points.mat','dam_points','refGK_LLcornerGLOBAL','dam_height_max','-v7.3')
end


%%% basins_selected
% parameters for the selected basin combination
if save_options.save_basins_selected == 1
    fields = {'dam_top_elev','dam_height','dam_heights','dam_elevations','dam_volumes',...
        'dam_axis_lengths','dam_axis_heights','basin_volumes','basin_areas',...
        'crit1_damvol_per_storage', 'crit2_area_per_storage', ...
        'crit3_share_wellsuited_cells', 'crit4_share_notsuited_cells'};
    for i_basins = 1:length(basins_selected)
        for i_field = 1:length(fields)
            basins_selected(i_basins).(fields{i_field}) = basins_selected(i_basins).(fields{i_field})(basins_selected(i_basins).dam_row_selected,:);
        end
    end
    fields = setdiff(fieldnames(basins_selected),...
        {'id', 'x_coord', 'y_coord', 'refGK_dam_coords', ...
        'dem','acc','dam_top_elev', 'dam_height','dam_axis_lengths',...
        'dam_axis_lengths_segments','dam_axis_heights',...
        'refGK_BASIN_x','refGK_BASIN_y',...
        'refBASIN_depth_dam_crosssection', 'refBASIN_depth_wall_crosssection',...
        'dam_heights','dam_elevations','dam_volumes','basin_volumes','basin_areas',...
        'crit1_damvol_per_storage', 'crit2_area_per_storage', ...
        'crit3_share_wellsuited_cells', 'crit4_share_notsuited_cells',...
        'ShA_dam', 'ShA_wall','curve_potential_dam_basin'});
    basins_selected = rmfield(basins_selected,fields);
    save('basins_selected.mat','basins_selected','refGK_LLcornerGLOBAL','-v7.3')
end


%%% save basins and dams as ascii-file
if save_options.save_basins_as_ascii == 1
    fun_9n1_save_basins_as_ascii    
end


%%% save basin characteristics as xlsx
if save_options.save_curves_as_excel == 1
    fun_9n2_save_basins_characteristics_as_xlsx    
end