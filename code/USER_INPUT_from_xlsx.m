function [save_options,plot_options,raster_selected] = USER_INPUT_from_xlsx(name_directory,name_file)
% ######  user input for basin location determination  ######
% summary of all user inputs which are required for the determination of
% retention basin locations (data is read from xlsx-file)
%
%
% Oinput variables and files:
%                     "xlsx_name" in "name_directory"
% functions: -
% called by m-codes: LOCASIN.m
%
%
% Author: Sonja Tesschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019


%% directories

global path_data path_result


% read excel-files
xlsx_name = 'copy_parameters.xlsx';
cd(name_directory)
copyfile(name_file,xlsx_name)
tbl_input_data = readtable(xlsx_name,'sheet','required_input_data');
tbl_exclude_dam = readtable(xlsx_name,'sheet','optional_exclude_dam');
tbl_exclude_basin = readtable(xlsx_name,'sheet','optional_exclude_basin');
tbl_parameters = readtable(xlsx_name,'sheet','parameters');
tbl_raster_selected = readtable(xlsx_name,'sheet','raster_selected');
delete(xlsx_name)


%directory with the input and output data
path_data = char(tbl_input_data{1,2});
path_result = char(tbl_input_data{2,2});
if exist(path_result,'dir')~=7
    mkdir(path_result);
end


%% input data

global grids_required info_exclude_dam info_exclude_basin

%%% REQUIRED %%%
% names of raster ascii-files which are included (located in path_data)
grids_required.dem = char(tbl_input_data{3,2});                % digital elevation model
grids_required.dem_fill = char(tbl_input_data{4,2});           % digital elevation model with filled sinks
grids_required.dir = char(tbl_input_data{5,2});                % flow directions
grids_required.acc = char(tbl_input_data{6,2});                % flow accumulation

grids_required.x_range =  [str2num(string(tbl_input_data{7,2})),str2num(string(tbl_input_data{8,2}))];
grids_required.y_range = [str2num(string(tbl_input_data{9,2})),str2num(string(tbl_input_data{10,2}))];


%%% OPTIONAL %%%
% info_exclude_dam: spatial description of the positions where a dam cannot
%                   be located (the basin is allowed on that position)
% info_exclude_basin: spatial description of the positions where no dam or
%                    basin is allowed (e.g. land use restrictions)

% exclude_dam
for i_data = 1:size(tbl_exclude_dam,1)
   field_name =  char(tbl_exclude_dam{i_data,1});
  if i_data==1 || isfield(info_exclude_dam,field_name)==0
       % name
       info_exclude_dam.(field_name).name = char(tbl_exclude_dam{i_data,2});
       % well-suited
       if isempty(tbl_exclude_dam{i_data,3})==0
           info_exclude_dam.(field_name).well_suited = str2num(string(tbl_exclude_dam{i_data,3}));
       end
       % not_suited
       if isempty(tbl_exclude_dam{i_data,4})==0
           info_exclude_dam.(field_name).not_suited = str2num(string(tbl_exclude_dam{i_data,4}));
       end
       % include
       if isempty(tbl_exclude_dam{i_data,5})==0
           info_exclude_dam.(field_name).include = str2num(string(tbl_exclude_dam{i_data,5}));
       end
       % exclude
       if isempty(tbl_exclude_dam{i_data,6})==0
           info_exclude_dam.(field_name).exclude{1} = str2num(string(tbl_exclude_dam{i_data,6}));
           info_exclude_dam.(field_name).exclude_threshold{1} = tbl_exclude_dam{i_data,7};
       end
   else
       % only the exclude-field can be defined more than once
       info_exclude_dam.(field_name).exclude{end+1} = str2num(string(tbl_exclude_dam{i_data,6}));
       info_exclude_dam.(field_name).exclude_threshold{end+1} = tbl_exclude_dam{i_data,7};
   end    
end

% exclude_basin
for i_data = 1:size(tbl_exclude_basin,1)
   field_name =  char(tbl_exclude_basin{i_data,1});
   if i_data==1 || isfield(info_exclude_basin,field_name)==0
       % name
       info_exclude_basin.(field_name).name = char(tbl_exclude_basin{i_data,2});
       % well-suited
       if isempty(tbl_exclude_basin{i_data,3})==0
           info_exclude_basin.(field_name).well_suited = str2num(string(tbl_exclude_basin{i_data,3}));
       end
       % not_suited
       if isempty(tbl_exclude_basin{i_data,4})==0
           info_exclude_basin.(field_name).not_suited = str2num(string(tbl_exclude_basin{i_data,4}));
       end
       % include
       if isempty(tbl_exclude_basin{i_data,5})==0
           info_exclude_basin.(field_name).include = str2num(string(tbl_exclude_basin{i_data,5}));
       end
       % exclude
       if isempty(tbl_exclude_basin{i_data,6})==0
           info_exclude_basin.(field_name).exclude{1} = str2num(string(tbl_exclude_basin{i_data,6}));
           info_exclude_basin.(field_name).exclude_threshold{1} = str2num(string(tbl_exclude_basin{i_data,7}));
       end
   else
       % only the exclude-field can be defined more than once
       info_exclude_basin.(field_name).exclude{end+1} = str2num(string(tbl_exclude_basin{i_data,6}));
       info_exclude_basin.(field_name).exclude_threshold{end+1} = str2num(string(tbl_exclude_basin{i_data,7}));
   end    
end


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
dam_height_max = tbl_parameters{1,2};
dam_height_min = tbl_parameters{2,2}; %0.2;%[3];%,3,2,1];  5
dam_height_buffer = tbl_parameters{3,2}; %.05;

% dam length in (m)
dam_length_max = tbl_parameters{4,2};   % the value must not be larger than the input grid size
exclude_longer_dams = tbl_parameters{5,2};   
    % 1: exclude dams, which are longer than dam_length_max
    % 0: do not exclude any dams because of its lengths
 
% basin volume in (m³)
basin_volume_max = tbl_parameters{6,2};
basin_volume_min = tbl_parameters{7,2};

% specific volume of the basin (only possible if flow accumulations were
% determined for the whole catchment)
% --> define empty if it should not be considered
sV_max = tbl_parameters{8,2};
sV_min = tbl_parameters{9,2};

% dam characteristics
dam_slope_m = tbl_parameters{10,2};         % slope in dam_slope_m per meter height
dam_crest_width = tbl_parameters{11,2};     % (m)


%%%%%   thresholds and discretization   %%%%%

% threshold for the definition of rivers from flow accumulation grid:
thresh = tbl_parameters{12,2};

% decide if all dam heights should be considered or only the dam axis
% heights (1: all dam points, 2: only dam axis)
limit_dam_height = tbl_parameters{13,2};

% evaluation of different dam heights
dam_dist_eval = tbl_parameters{14,2};   % step size between minimum and maximum dam height

% discretization level for curves (water depth - volume - area)
discretization_number = tbl_parameters{15,2};

% definition of the level of detail of the analysis
% (neighbors_exclude_distance = 0: every river point is analyzed,
% (neighbors_exclude_distance = 20: river points 20 m upstream and downstream of an
%  analyzed river points are excluded for further analysis)
neighbors_exclude_distance = tbl_parameters{16,2};


%%%%%   basin evaluation   %%%%%

% parameters of the objective function
w1_damVolume_per_basinVolume = tbl_parameters{17,2};
w2_basinArea_per_basinVolume = tbl_parameters{18,2};
w3_share_well_suited = tbl_parameters{19,2};
w4_share_not_suited = tbl_parameters{20,2};



%% debug mode and save options
global debug_on  save_memory

% definition, if debugging is on (plotting of figures)
debug_on = 0;

% decide whether to save memory (delete river_points, if exit_codes ~= 0)
save_memory = 1;

% definition, which output is saved
save_grids = tbl_parameters{21,2};
save_river_points = tbl_parameters{22,2};
save_dam_points = tbl_parameters{23,2};
save_basins_selected = tbl_parameters{24,2};
save_basins_as_ascii = tbl_parameters{25,2};
save_curves_as_excel = tbl_parameters{26,2};

save_options = struct('save_grids',save_grids,...
    'save_river_points',save_river_points,...
    'save_dam_points',save_dam_points,...
    'save_basins_selected',save_basins_selected,...
    'save_basins_as_ascii',save_basins_as_ascii,...
    'save_curves_as_excel',save_curves_as_excel);


% definition, which output is plotted
plot_exitcodes = tbl_parameters{27,2};
plot_spatial_overview = tbl_parameters{28,2};
plot_factsheet_p1 = tbl_parameters{29,2};
plot_factsheet_p2 = tbl_parameters{30,2};
plot_dam_comparison = tbl_parameters{31,2};
plot_curve_comparison = tbl_parameters{32,2};
plot_visibility = tbl_parameters{33,2};

plot_options = struct('plot_exitcodes',plot_exitcodes,...
    'plot_spatial_overview',plot_spatial_overview,...
    'plot_factsheet_p1',plot_factsheet_p1,...
    'plot_factsheet_p2',plot_factsheet_p2,...
    'plot_dam_comparison',plot_dam_comparison,...
    'plot_curve_comparison',plot_curve_comparison,...
    'plot_visibility',plot_visibility);


% define characteristic of raster for background map
raster_selected.name = char(tbl_raster_selected.name(1));
raster_selected.legend = tbl_raster_selected.legend';
raster_selected.color = [tbl_raster_selected.color_r,tbl_raster_selected.color_g,tbl_raster_selected.color_b];


