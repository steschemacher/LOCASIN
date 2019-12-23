function [save_options,plot_options,raster_selected] = USER_INPUT

% ######  user input for basin location determination  ######
% summary of all user inputs which are required for the determination of
% retention basin locations
%
% Other m-files required: -
% functions: -
% MAT-files required: -
% called by m-codes: LOCASIN.m
%
% Author: Sonja Tesschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019


%% directories
global path_data path_result

%directories of the input and output data
path_data = 'DIRECTORY IN WHICH THE DATA IS LOCATED';
path_result = 'DIRECTORY IN WHICH THE RESULTS ARE SAVED';
if exist(path_result,'dir')~=7
    mkdir(path_result);
end

%% input data
global grids_required info_exclude_dam info_exclude_basin

%%% REQUIRED %%%
% names of raster ascii-files which are included (located in path_data)
grids_required.dem = 'dem5.txt';                % digital elevation model
grids_required.dem_fill = 'dem_fill5.txt';      % digital elevation model with filled sinks
grids_required.dir = 'dir5.txt';                % flow directions
grids_required.acc = 'acc5.txt';                % flow accumulation

grids_required.x_range = [];
grids_required.y_range = [];


%%% OPTIONAL %%%
% info_exclude_dam: spatial description of the positions where a dam cannot
%                   be located (the basin is allowed on that position)
% info_exclude_basin: spatial description of the positions where no dam or
%                    basin is allowed (e.g. land use restrictions)
info_exclude_dam.river_buffer.name = 'buffer5.txt';     % real river net with buffer zone
info_exclude_dam.river_buffer.include = 70;             % include: numbers which are included, exclude: numbers which are excluded
info_exclude_basin.landuse.name = 'use5.txt';           % land use grid
info_exclude_basin.landuse.exclude{1} = [5,6];          % include: numbers which are included, exclude: numbers which are excluded
info_exclude_basin.landuse.exclude_threshold{1} = 100;  % threshold for exclusion cells in basins (in m²)
info_exclude_basin.landuse.exclude{2} = [4];            % include: numbers which are included, exclude: numbers which are excluded
info_exclude_basin.landuse.exclude_threshold{2} = 1000; % threshold for exclusion cells in basins (in m²)
info_exclude_basin.landuse.well_suited = [1];    % well_suited: numbers which are well suited for basins
info_exclude_basin.landuse.not_suited = [5,6];      % not_suited: numbers which are not suited for basins



%% parameters to define the dam and basin characteristics
% (thresholds, defintions, ...)
global thresh dam_height_max dam_height_min dam_length_max exclude_longer_dams
global basin_volume_max basin_volume_min dam_crest_width dam_slope_m
global dam_dist_eval discretization_number limit_dam_height
global neighbors_exclude_distance
global w1_damVolume_per_basinVolume w2_basinArea_per_basinVolume
global w3_share_well_suited w4_share_not_suited
global sV_min sV_max dam_height_buffer

%%%%%   limits for dams and basins   %%%%%

% maximum dam heigths which are analyzed:
dam_height_max = 6;
dam_height_min = 1;
dam_height_buffer = .05;

% dam length in (m)
dam_length_max = 70;   % the value must not be larger than the input grid size
exclude_longer_dams = 0;
    % 1: exclude dams, which are longer than dam_length_max
    % 0: do not exclude any dams because of its lengths

% basin volume in (m³)
basin_volume_max = 500000;
basin_volume_min = 5000;

% specific volume of the basin (only possible if flow accumulations were
% determined for the whole catchment)
% --> define empty if it should not be considered
sV_max = 30;
sV_min = 10;

% dam characteristics
dam_slope_m = 2;         % slope in dam_slope_m per meter height
dam_crest_width = 3;     % (m)


%%%%%   thresholds and discretization   %%%%%

% threshold for the definition of rivers from flow accumulation grid:
thresh = [];

% decide if all dam heights should be considered or only the dam axis
% heights (1: all dam points, 2: only dam axis)
limit_dam_height = 1;

% evaluation of different dam heights
dam_dist_eval = .2; % step size between minimum and maximum dam height

% discretization level for curves (water depth - volume - area)
discretization_number = 41;

% definition of the level of detail of the analysis
% (neighbors_exclude_distance = 0: every river point is analyzed,
% (neighbors_exclude_distance = 20: river points 20 m upstream and downstream of an
%  analyzed river points are excluded for further analysis)
neighbors_exclude_distance = 0;


%%%%%   basin evaluation   %%%%%

% parameters of the objective function
w1_damVolume_per_basinVolume = .3;
w2_basinArea_per_basinVolume = .3;
w3_share_well_suited = .2;
w4_share_not_suited = .2;



%% debug mode and save options
global debug_on  save_memory

% definition, if debugging is on (plotting of figures)
debug_on = 0;

% decide whether to save memory (delete river_points, if exit_codes ~= 0)
save_memory = 1;

% definition, which output is saved
save_grids = 1;
save_river_points = 1;
save_dam_points = 1;
save_basins_selected = 1;
save_basins_as_ascii = 1;
save_curves_as_excel = 1;

save_options = struct('save_grids',save_grids,...
    'save_river_points',save_river_points,...
    'save_dam_points',save_dam_points,...
    'save_basins_selected',save_basins_selected,...
    'save_basins_as_ascii',save_basins_as_ascii,...
    'save_curves_as_excel',save_curves_as_excel);

% definition, which output is plotted
plot_exitcodes = 1;
plot_spatial_overview = 1;
plot_factsheet_p1 = 1;
plot_factsheet_p2 = 1;
plot_dam_comparison = 1;
plot_curve_comparison = 1;
plot_visibility = 0;   % 1: figures are visible, 0: figures are invisible

plot_options = struct('plot_exitcodes',plot_exitcodes,...
    'plot_spatial_overview',plot_spatial_overview,...
    'plot_factsheet_p1',plot_factsheet_p1,...
    'plot_factsheet_p2',plot_factsheet_p2,...
    'plot_dam_comparison',plot_dam_comparison,...
    'plot_curve_comparison',plot_curve_comparison,...
    'plot_visibility',plot_visibility);


% define characteristic of raster for background map
raster_selected.name = 'landuse';
raster_selected.legend = {'pasture','cropland','forest','track','road','sealed'};
raster_selected.color = [[86, 129, 45]/255;
    [144, 131, 12]/255;
    [43, 64, 27]/255;
    [102, 105, 99]/255;
    [177, 108, 37]/255;
    [141, 36, 37]/255];


