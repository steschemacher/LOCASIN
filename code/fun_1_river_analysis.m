function fun_1_river_analysis
% ######  river analysis  ######
% function to analyze and characterize river points
%
% functions:    fun_1n1_ASCIimport.m
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019


global grids_required path_data thresh info_exclude_dam info_exclude_basin dam_length_max 
global river_points grids grids_exclude point_boundaries dam_height_max dam_slope_m
global refGK_LLcornerGLOBAL
global basin_volume_max basin_volume_min sV_min sV_max
global river_longitudinal_section

%% loading input grids (output ascii-files from arcmap)

fields = fieldnames(grids_required);
fields = setdiff(fields,{'x_range','y_range'});

for i_field = 1:length(fields)
    
    % load data
    if isfield(grids_required,'x_range')==1 && isempty(grids_required.x_range)==0
        tmp = fun_1n1_ASCIimport(path_data,grids_required.(fields{i_field}),grids_required.x_range,grids_required.y_range);
    else
        tmp = fun_1n1_ASCIimport(path_data,grids_required.(fields{i_field}),[],[]);
    end
    grids.(fields{i_field}) = tmp.data;
    
    % include raster characteristics
    if i_field == 1
        grids.ncols = tmp.ncols;
        grids.nrows = tmp.nrows;
        grids.xll = tmp.xll;
        grids.yll = tmp.yll;
        grids.cellsize = tmp.cellsize;
        grids.x = tmp.x;
        grids.y = tmp.y;
        grids.X = tmp.X;
        grids.Y = tmp.Y;
    end
    clear tmp
end

% define threshold (if necessary)
if (isempty(thresh)==1 || isnan(thresh)==1) && isempty(sV_max)==0
    thresh = basin_volume_min/sV_max/(grids.cellsize)^2*1e6/1e3 / 5;    
end

% determine river network
RIV = grids.acc;
RIV(RIV==-9999) = 0;                    % set nodata-values to 0
RIV(RIV<=thresh) = 0;                   % set all rivers to 1
RIV(RIV>thresh) = 1;                    % set all rivers to 1
grids.riv = RIV;
clear RIV

% definition of the required boundary for a river point for the analysis
point_boundaries = round(dam_length_max/2/grids.cellsize + 1.5*max(dam_height_max)*dam_slope_m/grids.cellsize);

% reference point of the lower left corner in Gauss-Krueger-coordinates
refGK_LLcornerGLOBAL = [grids.xll;grids.yll];


%% analyze the river net

% rows and cols of river points
[row_coord,col_coord] = size(grids.riv);
COL_nums = repmat(1:col_coord,row_coord,1);       % matric with x-coordinates
ROW_nums = repmat((1:row_coord)',1,col_coord);    % matrix with y-coordinates
% table with river points and coordinates
clear river_points_tmp
river_points_tmp(:,2:4) = [reshape(COL_nums,[],1) reshape(ROW_nums,[],1) reshape(grids.riv,[],1)];
river_points_tmp = river_points_tmp(river_points_tmp(:,4)==1,1:3);
river_points_tmp(:,1) = 1:size(river_points_tmp,1);

% indices of the river points
index_river = sub2ind(size(grids.riv),river_points_tmp(:,3),river_points_tmp(:,2));


%%% Nachfolger / downstream neighbor
id_grid = nan(size(grids.riv));
id_grid(index_river) = ind2sub(size(grids.riv),river_points_tmp(:,1));% (1:numel(grids.riv),size(grids.riv,1),size(grids.riv,2));
nachfolger_grid = nan(size(grids.riv));

% rechts-oben
id_grid_move = nan(size(grids.riv));
id_grid_move(2:end-1,2:end-1) = id_grid(1:end-2,3:end);
nachfolger_grid(grids.dir==128) = id_grid_move(grids.dir==128);
% rechts
id_grid_move = nan(size(grids.riv));
id_grid_move(2:end-1,2:end-1) = id_grid(2:end-1,3:end);
nachfolger_grid(grids.dir==1) = id_grid_move(grids.dir==1);
% rechts-unten
id_grid_move = nan(size(grids.riv));
id_grid_move(2:end-1,2:end-1) = id_grid(3:end,3:end);
nachfolger_grid(grids.dir==2) = id_grid_move(grids.dir==2);
% unten
id_grid_move = nan(size(grids.riv));
id_grid_move(2:end-1,2:end-1) = id_grid(3:end,2:end-1);
nachfolger_grid(grids.dir==4) = id_grid_move(grids.dir==4);
% links-unten
id_grid_move = nan(size(grids.riv));
id_grid_move(2:end-1,2:end-1) = id_grid(3:end,1:end-2);
nachfolger_grid(grids.dir==8) = id_grid_move(grids.dir==8);
% links
id_grid_move = nan(size(grids.riv));
id_grid_move(2:end-1,2:end-1) = id_grid(2:end-1,1:end-2);
nachfolger_grid(grids.dir==16) = id_grid_move(grids.dir==16);
% links-oben
id_grid_move = nan(size(grids.riv));
id_grid_move(2:end-1,2:end-1) = id_grid(1:end-2,1:end-2);
nachfolger_grid(grids.dir==32) = id_grid_move(grids.dir==32);
% oben
id_grid_move = nan(size(grids.riv));
id_grid_move(2:end-1,2:end-1) = id_grid(1:end-2,2:end-1);
nachfolger_grid(grids.dir==64) = id_grid_move(grids.dir==64);


%%% Vorgaenger / upstream neighbor

nachfolger_array = [river_points_tmp(:,1),nachfolger_grid(index_river)];
clear vorgaenger_array vorgaenger_cell
[C,ia,~] = unique(nachfolger_array(:,2));
vorgaenger_array(1:size(nachfolger_array,1),1) = 0;
vorgaenger_array(C(isnan(C)==0),1) = nachfolger_array(ia(isnan(C)==0),1);
missing_vorgaenger = nachfolger_array(ismember(nachfolger_array(:,1),setdiff(nachfolger_array(isnan(nachfolger_array(:,2))==0,1),vorgaenger_array)),:);
while isempty(missing_vorgaenger)==0
    [C,ia,~] = unique(missing_vorgaenger(:,2));
    vorgaenger_array(:,end+1) = 0;
    vorgaenger_array(C(isnan(C)==0),end) = missing_vorgaenger(ia(isnan(C)==0),1);
    missing_vorgaenger = missing_vorgaenger(ismember(missing_vorgaenger(:,1),setdiff(missing_vorgaenger(isnan(missing_vorgaenger(:,2))==0,1),missing_vorgaenger)),:);
end

vorgaenger_cell = num2cell(vorgaenger_array,2);
vorgaenger_cell = cellfun(@(x)x(x ~= 0), vorgaenger_cell, 'UniformOutput', false);


id_river_points(:,1) = sub2ind(size(grids.dem),river_points_tmp(:,3),river_points_tmp(:,2));


% generate dummy field dem_fill, if it does not exist (which is later
% deleted again)
if isfield(grids,'dem_fill')==0 && isempty(grids.dem_fill)==0
    grids.dem_fill = grids.dem;
    delete_dem_fill = 1;
else
    delete_dem_fill = 0;
end

% struct with river points and characteristics
river_points = struct('id',num2cell(river_points_tmp(:,1)),...
    'id_grid',num2cell(id_river_points(:,1)),...
    'row',num2cell(river_points_tmp(:,3)),...
    'col',num2cell(river_points_tmp(:,2)),...
    'x_coord',num2cell(grids.X(index_river)),...
    'y_coord',num2cell(grids.Y(index_river)),...
    'dem',num2cell(grids.dem(index_river)),...   
    'dem_fill',num2cell(grids.dem_fill(index_river)),... 
    'acc',num2cell(grids.acc(index_river)),...
    'vorgaenger',vorgaenger_cell,...
    'nachfolger',num2cell(nachfolger_grid(index_river)),...
    'alle_vorgaenger',[],...
    'exit_code',0);

clear river_points_tmp


%% determine all upstream neighbors

start_points_all =  find(arrayfun(@(river_points) isempty(river_points.vorgaenger),river_points));
end_points_all = find(arrayfun(@(river_points) isnan(river_points.nachfolger),river_points));
start_points = setdiff(start_points_all,end_points_all);
end_points = setdiff(end_points_all,start_points_all);

clear start_points_all end_points_all

for i_start = 1:length(start_points)
    river_longitudinal_section{i_start} = [];
    i_location_upstream = start_points(i_start);    
    i_location = river_points(i_location_upstream).nachfolger;
    while sum(ismember(end_points,i_location))~=1  
       river_longitudinal_section{i_start} =  [river_longitudinal_section{i_start},i_location_upstream];
       i_location_upstream = river_points(i_location).id;       
       i_location = river_points(i_location).nachfolger;
    end  
end


%% exclude cells with the wrong specific volume

if isempty(sV_min)==0
    points_sV_basin_volume_max = basin_volume_max./([river_points.acc]*(grids.cellsize)^2/1e6)/1e3;
    points_exclude = find(points_sV_basin_volume_max < sV_min);
    if isempty(points_exclude)==0
        points_excluded = find([river_points.exit_code]>0);
        points_exclude = setdiff(points_exclude,points_excluded);
        [river_points(points_exclude).exit_code] = deal(3);
    end
end
if isempty(sV_max)==0
    points_sV_basin_volume_min = basin_volume_min./([river_points.acc]*(grids.cellsize)^2/1e6)/1e3;
    points_exclude = find(points_sV_basin_volume_min > sV_max);
    if isempty(points_exclude)==0
        points_excluded = find([river_points.exit_code]>0);
        points_exclude = setdiff(points_exclude,points_excluded);
        [river_points(points_exclude).exit_code] = deal(4);
    end
end


%% exclude river points due to spatial data --> exit_code = 1

% raster index of the river points (raster sizes need to be equal!!!)
index_river = sub2ind(size(grids.riv),[river_points.row],[river_points.col]);
info_exclude_dam.borders.exclude{1} = 1;   % guarantee that the dams fit into the raster

% load data as basis for the exclusion of river points (info_exclude_dam)
fields = fieldnames(info_exclude_dam);

for i_field = 1:length(fields)
    
    if i_field < length(fields)   % for external data
    % load data
    if isfield(grids_required,'x_range')==1 && isempty(grids_required.x_range)==0
        tmp = fun_1n1_ASCIimport(path_data,info_exclude_dam.(fields{i_field}).name,grids_required.x_range,grids_required.y_range);
    else
        tmp = fun_1n1_ASCIimport(path_data,info_exclude_dam.(fields{i_field}).name,[],[]);
    end
    grids_exclude.(fields{i_field}) = tmp.data;
    
    % include raster characteristics
    if i_field == 1
        grids_exclude.ncols = tmp.ncols;
        grids_exclude.nrows = tmp.nrows;
        grids_exclude.xll = tmp.xll;
        grids_exclude.yll = tmp.yll;
        grids_exclude.cellsize = tmp.cellsize;
        grids_exclude.x = tmp.x;
        grids_exclude.y = tmp.y;
        grids_exclude.X = tmp.X;
        grids_exclude.Y = tmp.Y;
    end
    clear tmp
    else    % exlcude data which is too close to the borders of the rasters
        grids_exclude.borders = zeros(size(grids_exclude.(fields{i_field-1})));
        grids_exclude.borders(1:point_boundaries,:) = 1;
        grids_exclude.borders(end-point_boundaries:end,:) = 1;
        grids_exclude.borders(:,1:point_boundaries) = 1;
        grids_exclude.borders(:,end-point_boundaries:end) = 1;
    end
    
    % select points which have to be excluded (or negative of included)
    if isfield(info_exclude_dam.(fields{i_field}),'exclude')==1 && isempty(info_exclude_dam.(fields{i_field}).exclude{1})==0
        grid_tmp = grids_exclude.(fields{i_field});
        for i_exclude = 1:length(info_exclude_dam.(fields{i_field}).exclude)
            lines_exclude = ismember(grid_tmp(index_river),info_exclude_dam.(fields{i_field}).exclude{i_exclude});
            lines_exclude = [river_points(lines_exclude==1).id];
            lines_excluded = find([river_points.exit_code]>0);
            lines_exclude = setdiff(lines_exclude,lines_excluded);
            [river_points(lines_exclude).exit_code] = deal(1);
        end
    elseif isfield(info_exclude_dam.(fields{i_field}),'include')==1 && isempty(info_exclude_dam.(fields{i_field}).include)==0
        grid_tmp = grids_exclude.(fields{i_field});
        lines_include = ismember(grid_tmp(index_river),info_exclude_dam.(fields{i_field}).include);
        lines_include = [river_points(lines_include~=1).id];
        [river_points(lines_include).exit_code] = deal(1);
    end   
    
end

% load data as basis for the exclusion of river points (info_exclude_basin)
fields = fieldnames(info_exclude_basin);

for i_field = 1:length(fields)
    % load data
    if isfield(grids_required,'x_range')==1 && isempty(grids_required.x_range)==0
        tmp = fun_1n1_ASCIimport(path_data,info_exclude_basin.(fields{i_field}).name,grids_required.x_range,grids_required.y_range);
    else
        tmp = fun_1n1_ASCIimport(path_data,info_exclude_basin.(fields{i_field}).name,[],[]);
    end
    grids_exclude.(fields{i_field}) = tmp.data;
    
    % include raster characteristics
    if i_field == 1
        grids_exclude.ncols = tmp.ncols;
        grids_exclude.nrows = tmp.nrows;
        grids_exclude.xll = tmp.xll;
        grids_exclude.yll = tmp.yll;
        grids_exclude.cellsize = tmp.cellsize;
        grids_exclude.x = tmp.x;
        grids_exclude.y = tmp.y;
        grids_exclude.X = tmp.X;
        grids_exclude.Y = tmp.Y;
    end
    clear tmp
    
    
    % select points which have to be excluded (or negative of included)
    if isfield(info_exclude_basin.(fields{i_field}),'exclude')==1 && isempty(info_exclude_basin.(fields{i_field}).exclude)==0
        grid_tmp = grids_exclude.(fields{i_field});
        for i_exclude = 1:length(info_exclude_basin.(fields{i_field}).exclude)
            lines_exclude = ismember(grid_tmp(index_river),info_exclude_basin.(fields{i_field}).exclude{i_exclude});
            lines_exclude = [river_points(lines_exclude==1).id];
            lines_excluded = find([river_points.exit_code]>0);
            lines_exclude = setdiff(lines_exclude,lines_excluded);
            [river_points(lines_exclude).exit_code] = deal(2);
        end
    elseif isfield(info_exclude_basin.(fields{i_field}),'include')==1 && isempty(info_exclude_basin.(fields{i_field}).include)==0
        grid_tmp = grids_exclude.(fields{i_field});
        lines_include = ismember(grid_tmp(index_river),info_exclude_basin.(fields{i_field}).include);
        lines_include = [river_points(lines_include~=1).id];
        lines_included = find([river_points.exit_code]>0);
        lines_include = setdiff(lines_include,lines_included);
        [river_points(lines_include).exit_code] = deal(2);
    end   
    
end


%% delete fields which are not needed for further analyses
% (to decrease memory usage)

if delete_dem_fill ==1
    river_points = rmfield(river_points,'dem_fill');
end
% river_points = rmfield(river_points,'dir');
grids = rmfield(grids,'acc');
grids = rmfield(grids,'dem_fill');
grids = rmfield(grids,'dir');
grids = rmfield(grids,'riv');
grids = rmfield(grids,'X');
grids = rmfield(grids,'Y');

