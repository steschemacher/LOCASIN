% ######  tool for basin location determination  ######
% analysis of the topography in order to find basin locations
%
%
% input variables and files:
%               define_input_directory_and_file.txt
%               (line 1: directory, line 2: file)
%               read_from_xlsx (1/0: yes/no)
%
% functions:    USER_INPUT.m / USER_INPUT_from_xlsx.m (user_input.xlsx)
%               fun_1_river_analysis.m
%                   (fun_ASCIIimport.m)
%               fun_2_determine_shortest_dam.m
%               fun_3_determine_river_points_in_basin.m
%               fun_4_determine_dam_characteristics.m
%               fun_5_determine_basin_area_and_volume.m
%               fun_6_evaluate_basins.m
%               fun_7_selection_of_basin_combinations
%               fun_8_determine_depth_storage_area_curves.m
%               fun_9_summary_and_save_output.m
%               fun_10_plot_results.m
%
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019


clear all
close all
clc

% decision whether to read data from xlsx or m-file
read_from_xlsx = 1;

global river_points  river_points_all save_memory basins_selected dam_points


if read_from_xlsx==1
    info_input = textread('define_input_directory_and_file.txt','%s');
    name_directory = info_input{1};
    try
        cd(name_directory)
    catch
        error(sprintf('The defined input directory is not existing:\n %s',name_directory))
    end
end
diary('logfile.txt')

try
    %-- define status figure ---------------------------------------------%
    logo = imread('locasin.jpg');
    overview = figure('name','status','NumberTitle','off','visible','off');
    s = subplot(1,1,1);
    imshow(logo)
    size_position = get(overview,'Position');
    set(overview,'Position',[size_position(1),size_position(2),min(size(logo,2),round(size_position(3)*.8)),min(size(logo,2),round(size_position(3)*.8))/size(logo,2)*size(logo,1)]);
    set(s,'Position',[0,0,1,1])
    overview_fontsize = (min(size(logo,2),round(size_position(3)*.8))/size(logo,2)*size(logo,1))/30;
    set(overview,'visible','on')
    pause(1)
    %---------------------------------------------------------------------%
    
    
    %% PREPROCESSING
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.45,.5,.05],'String','analyzing river network...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start_total = now();
    %---------------------------------------------------------------------%
    
    
    %#####################################################################%
    %%%   read user input                                               %%%
    %#####################################################################%
    % --> selection of the desired input-data-format
    
    if read_from_xlsx ==1
        name_file = info_input{2};
        % read user input from user_input.xlsx-file
        [save_options,plot_options,raster_selected] = USER_INPUT_from_xlsx(name_directory,name_file);
    else
        % read user input from USER_INPUT.m-function
        [save_options,plot_options,raster_selected] = USER_INPUT;
    end
    
    
    
    
    %#####################################################################%
    %%%   analyze the river network                                     %%%
    %#####################################################################%
    
    fun_1_river_analysis;
    
    
    % copy river_points-struct and thin out the variables
    river_points_all = river_points;
    fields = setdiff(fieldnames(river_points_all),{'id','id_grid','x_coord','y_coord',...
        'exit_code','row','col','dem','dem_fill','acc','vorgaenger','nachfolger'});
    river_points_all = rmfield(river_points_all,fields);
    
    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start_total)*24*60;
    fprintf('Preprocessing completed: %3.3f min\n',duration);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.45,.5,.05],'String',sprintf('finished: %7.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %---------------------------------------------------------------------%
    
    
    
    %% DAM ANALYSIS
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.4,.5,.05],'String','analyzing dam orientations...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start = now();
    %---------------------------------------------------------------------%
    
    
    % select points which have not been excluded before
    if save_memory == 0
        points_to_analyze = [river_points([river_points.exit_code]<=0).id];
    else
        river_points = river_points([river_points.exit_code]<=0);
        points_to_analyze = 1:length(river_points);
    end
    
    for i_point = points_to_analyze        
        if river_points(i_point).exit_code <= 0
            
            i_point_all = find([river_points_all.id]==river_points(i_point).id);
            
            %#############################################################%
            %%%   detect dam endpoints (if possible)                    %%%
            %#############################################################%
            % determine the shortest dam at a river point to close the valley
            
            fun_2_determine_shortest_dam(i_point,i_point_all)
            
            
        end
    end
    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start)*24*60;
    fprintf('Dam analysis completed: %3.3f min\n',duration);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.4,.5,.05],'String',sprintf('finished: %7.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %---------------------------------------------------------------------%
    
    
    %% BASIN ANALYSIS
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.35,.5,.05],'String','analyzing basin locations...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start = now();
    %---------------------------------------------------------------------%
    
    
    % select points which have not been excluded before
    if save_memory == 0
        points_to_analyze = [river_points([river_points.exit_code]<=0).id];
    else
        river_points = river_points([river_points.exit_code]<=0);
        points_to_analyze = [river_points.id];
    end
    
    for id_analyze = points_to_analyze
        
        i_point = find([river_points.id]==id_analyze);
        i_point_all = find([river_points_all.id]==id_analyze);
        
        %#################################################################%
        %%%   determine river points within the potential basin         %%%
        %#################################################################%
        % determine all river points, which are potentially in the basin by
        % comparison of the elevation of the upstream river points with the
        % top elevation of the dam
        
        fun_3_determine_river_points_in_basin(i_point,i_point_all);
        
        if river_points(i_point).exit_code>0
            if save_memory == 1
                % exclude site
                river_points = river_points([river_points.exit_code]<=0);
            end
            % no further analysis required (no basin end point existent)
            continue
        end
        
        
        %#################################################################%
        %%%   determine dam characteristics and volume                  %%%
        %#################################################################%
        % determine dam characteristics for different heights: dam heights,
        % dam elevations, dam volumes, dam axis lengths, dam axis heights
        
        fun_4_determine_dam_characteristics(i_point);
        
        if river_points(i_point).exit_code>0
            if save_memory == 1
                % exclude site
                river_points = river_points([river_points.exit_code]<=0);
            end
            % no further analysis required (all pot. dams are too high)
            continue
        end
        
        
        %#################################################################%
        %%%   determine basin area and volume                           %%%
        %#################################################################%
        
        fun_5_determine_basin_area_and_volume(i_point,i_point_all);
        
        
        if river_points(i_point).exit_code>0
            if save_memory == 1
                % exclude site
                river_points = river_points([river_points.exit_code]<=0);
            end
            % no further analysis required (no basin possible)
            continue
        end
        
        
        %#################################################################%
        %%%   rate the determined basins                                %%%
        %#################################################################%
        
        fun_6_evaluate_basins(i_point,i_point_all);
        
        if river_points(i_point).exit_code>0
            if save_memory == 1
                % exclude site
                river_points = river_points([river_points.exit_code]<=0);
            end
            % no further analysis required (no basin possible)
            continue
        end
        
    end
    
    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start)*24*60;
    fprintf('Basin analysis completed: %3.3f min\n',duration);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.35,.5,.05],'String',sprintf('finished: %7.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %---------------------------------------------------------------------%
    
    
    %% BASIN COMBINATION
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.3,.5,.05],'String','selecting basin combination...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start = now();
    %---------------------------------------------------------------------%
    
    
    if save_memory == 0
        % select points which have not been excluded before
        points_to_analyze = [river_points([river_points.exit_code]<=0).id];
        dam_points = river_points(points_to_analyze);
    else
        dam_points = river_points;
    end
    
    
    %#####################################################################%
    %%%   selection of the best basin combination                       %%%
    %#####################################################################%
    
    if isempty(dam_points)==0
        basins_selected = fun_7_selection_of_basin_combination(dam_points);
    else
        disp('No basin site available.')
        % exclude saving-options with empty files
        save_dam_points = 0;
        save_basins_selected = 0;
        save_basins_as_ascii = 0;
        save_curves_as_excel = 0;
        % exclude plotting-options with empty plots
        plot_spatial_overview = 0;
        plot_factsheet_p1 = 0;
        plot_factsheet_p2 = 0;
        plot_dam_comparison = 0;
        plot_curve_comparison = 0;
    end
    
    
    %#####################################################################%
    %%%   determine curves for water depth, basin area and basin storage %%
    %#####################################################################%
    
    for i_point = 1:length(basins_selected)
        fun_8_determine_depth_storage_area_curves(i_point);
    end
    
    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start)*24*60;
    fprintf('Basin combination selected: %3.3f min\n',duration);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.3,.5,.05],'String',sprintf('finished: %7.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %---------------------------------------------------------------------%
    
    
    %% POSTPROCESSING
    
    %#####################################################################%
    %%%   summary of the results: save output-files                     %%%
    %#####################################################################%
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.25,.5,.05],'String','saving results...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start = now();
    %---------------------------------------------------------------------%
    
    
    %%% save results %%%
    fun_9_summary_and_save_output(save_options)
    
    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start)*24*60;
    fprintf('Results saved: %3.3f min\n',duration);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.25,.5,.05],'String',sprintf('finished: %7.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %---------------------------------------------------------------------%
    
    
    %#####################################################################%
    %%%   summary of the results: save plots                            %%%
    %#####################################################################%
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.2,.5,.05],'String','plotting results...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start = now();
    %---------------------------------------------------------------------%
    
    
    %%% plot figures %%%
    fun_10_plot_results(plot_options,raster_selected);
    
    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start)*24*60;
    duration_total = (time_end-time_start_total)*24*60;
    fprintf('Figures saved: %3.3f min\n',duration);
    fprintf('Total duration: %3.3f min\n',duration_total);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.2,.5,.05],'String',sprintf('finished: %7.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
        pause(1)
        annotation('textbox',[.45,.12,.5,.05],'String',sprintf('Total duration: %7.3f min\n',duration_total),'linestyle','none','fontsize',overview_fontsize)
    end
    %---------------------------------------------------------------------%
    diary off
    
catch ME
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.45,.12,.4,.05],'String','ERROR. Read logfile for more details.','linestyle','none','fontsize',overview_fontsize)
    end
    %---------------------------------------------------------------------%
    
    rethrow(ME)
    diary off
end

%% exit codes
%
%
% LOCASIN.m
%   -   no exclusions
%
% fun_1_river_analysis.m:
%   1:  river point too close to border or excluded by spatial data for dam
%       restrictions (e.g. buffer of river network) which was defined by
%       the user in USER_INPUT.m (name of grid and numbers to include or
%       exclude)
%   2:  exclusion by spatial data for basin restrictions (e.g. land use)
%       which was defined by the user in USER_INPUT.m (name of grid and
%       numbers to include or exclude)
%   3:  target basin volume is too small for target specific volume
%       --> (too) large rivers are excluded
%       [the maximum basin volume (user definition) would result in a
%       specific volume at the site, which is lower than the minimum
%       specific volume (user definition)]
%   4:  target basin volume is too large for target specific volume
%       --> (too) small rivers are exclued
%       [the minimum basin volume (user definition) would result in a
%       specific volume at the site, which is higher than the maximum
%       specific volume (user definition)]
%
% fun_2_determine_shortest_dam.m
%   5:  dam is excluded for further analysis, because it is a neighbor of
%       an already determined dam (user definition to reduce the
%       computational effort)
%   6:  there is no dam possible to close the valley
%       (for >= minimum height and length)
%   -1: dam height had to be reduced to close the valley
%
% fun_3_determine_river_points_in_basin.m
%   7:  no basin end point exists (river point)
%   -2: dam height had to be reduced, because the basin endpoint could not
%       be found within the extent of the DEM
%
% fun_4_determine_dam_characteristics.m
%   -   no exclusions
%
% fun_5_determine_basin_area_and_volume.m
%   8:  no upstream river point availabe for the basin area analysis
%   9:  only river points on the shoreline or dam could be found
%   10: basin touches the borders of the selected DEM
%       (usually occurs, if the dam does not close the valley totally, e.g.
%       if there is a conjunction within the basin)
%   -3: dam height had to be reduced, because the basin touched the borders
%       of the DEM
%
% fun_6_evaluate_basins.m
%   11: dam (points) too high
%   12: maximum dam height not high enough
%   13: dam axis is too long
%   14: dam volumes not large enough
%   15: dam volumes too large
%   16: too many exclusion-cells in the basin
%       (definition of the specific number by user)
%   17: specific volume too small
%   18: specific volume too large
%
% fun_7_selection_of_basin_combination.m
%   -   no exclusions
%
% fun_8_determine_depth_storage_area_curves.m
%   -   no exclusions
%
%
