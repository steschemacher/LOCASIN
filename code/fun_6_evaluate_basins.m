function fun_6_evaluate_basins(i_point,i_point_all)
% ######  evaluate dam heights  ######
% function to evaluate and exclude dam heights based on the dam and basin
% characteristics
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019


global river_points river_points_all grids_exclude info_exclude_basin grids
global dam_height_max dam_height_min basin_volume_max basin_volume_min
global sV_min sV_max river_longitudinal_section
global dam_length_max exclude_longer_dams dam_height_buffer
global refGK_LLcornerGLOBAL
global w1_damVolume_per_basinVolume w2_basinArea_per_basinVolume
global w3_share_well_suited w4_share_not_suited limit_dam_height

dam_heights_potential = river_points(i_point).dam_heights_eval;


%% exclude potential dams, which do not fit the criteria

if limit_dam_height == 1
    % limit of the maximum dam height should be applied for the whole dam
    dam_height_max_real = max(river_points(i_point).dam_heights,[],2);
    [dam_ids_exclude,~] = find(dam_height_max_real > max(dam_height_max) + dam_height_buffer);
elseif limit_dam_height == 2
    % limit of the maximum dam height should only be applied for the dam axis
    dam_height_max_real = max(river_points(i_point).dam_axis_heights,[],2);
    [dam_ids_exclude,~] = find(dam_height_max_real > max(dam_height_max) + dam_height_buffer);
end

dam_heights_potential(1,dam_ids_exclude) = NaN;
clear dam_ids_exclude

if sum(isnan(dam_heights_potential)==1) == length(dam_heights_potential)
    river_points(i_point).exit_code = 11; % dam points too high
    river_points_all(i_point_all).exit_code = 11; % dam points too high
    return
end

% minimum dam height
[dam_ids_exclude,~] = find(max(river_points(i_point).dam_axis_heights,[],2) < dam_height_min);
dam_heights_potential(1,dam_ids_exclude) = NaN;
clear dam_ids_exclude

if sum(isnan(dam_heights_potential)==1) == length(dam_heights_potential)
    river_points(i_point).exit_code = 12; % maximum dam height not high enough
    river_points_all(i_point_all).exit_code = 12; % maximum dam height not high enough
    return
end

% maximum dam length
if exclude_longer_dams == 1
    [dam_ids_exclude,~] = find(max(river_points(i_point).dam_axis_lengths,[],2) > dam_length_max);
    dam_heights_potential(1,dam_ids_exclude) = NaN;
    clear dam_ids_exclude
    
    if sum(isnan(dam_heights_potential)==1) == length(dam_heights_potential)
        river_points(i_point).exit_code = 13; % dam is too long
        river_points_all(i_point_all).exit_code = 13; % dam is too long
        return
    end
end

% minimum storage volume
[dam_ids_exclude,~] = find(river_points(i_point).basin_volumes(:,1) < basin_volume_min);
dam_heights_potential(1,dam_ids_exclude) = NaN;
clear dam_ids_exclude

if sum(isnan(dam_heights_potential)==1) == length(dam_heights_potential)
    river_points(i_point).exit_code = 14; % dam volumes not large enough
    river_points_all(i_point_all).exit_code = 14; % dam volumes not large enough
    return
end

% maximum storage volume
[dam_ids_exclude,~] = find(river_points(i_point).basin_volumes(:,1) > basin_volume_max);
dam_heights_potential(1,dam_ids_exclude) = NaN;
clear dam_ids_exclude

if sum(isnan(dam_heights_potential)==1) == length(dam_heights_potential)
    river_points(i_point).exit_code = 15; % dam volumes too large
    river_points_all(i_point_all).exit_code = 15; % dam volumes too large
end



%% exclude potential dams, which have too many restricted cells

% ids of refBASIN-raster
refBASIN = river_points(i_point).refBASIN_depth_dam_crosssection(:,:,1);
RASTER_ID = reshape(1:numel(refBASIN),size(refBASIN));


fields = fieldnames(info_exclude_basin);

ids_EXCLUDE{1:length(fields)} = RASTER_ID;
ids_WELL_SUITED = RASTER_ID;
ids_NOT_SUITED = RASTER_ID;


for i_field = 1:length(fields)
    
    % clip field-raster to the extent of the basin raster (refBASIN)
    grid_tmp = grids_exclude.(fields{i_field});
    grid_tmp = grid_tmp(river_points(i_point).refBASIN_LLcorner(2):...
        river_points(i_point).refBASIN_LLcorner(2)+size(refBASIN,1)-1,...
        river_points(i_point).refBASIN_LLcorner(1):...
        river_points(i_point).refBASIN_LLcorner(1)+size(refBASIN,2)-1);
    
    % select points which have to be excluded (or negative of included)
    if isfield(info_exclude_basin.(fields{i_field}),'exclude')==1 && isempty(info_exclude_basin.(fields{i_field}).exclude)==0
        for i_exclude = 1:length(info_exclude_basin.(fields{i_field}).exclude)
            lines_exclude = ismember(grid_tmp,info_exclude_basin.(fields{i_field}).exclude{i_exclude});        
            ids_EXCLUDE{i_field,i_exclude}(lines_exclude~=1) = NaN;
        end
    else 
        ids_EXCLUDE{i_field,1}(:,:) = NaN;
    end
    
    % select points which are WELL suited
    if isfield(info_exclude_basin.(fields{i_field}),'well_suited')==1 && isempty(info_exclude_basin.(fields{i_field}).well_suited)==0
        lines_wellsuited = ismember(grid_tmp,info_exclude_basin.(fields{i_field}).well_suited);
        ids_WELL_SUITED(lines_wellsuited~=1) = NaN;
    end
    
    % select points which are NOT suited
    if isfield(info_exclude_basin.(fields{i_field}),'not_suited')==1 && isempty(info_exclude_basin.(fields{i_field}).not_suited)==0
        lines_notsuited = ismember(grid_tmp,info_exclude_basin.(fields{i_field}).not_suited);
        ids_NOT_SUITED(lines_notsuited~=1) = NaN;
    end
    
    for i_exclude = 1:size(ids_EXCLUDE,2)
        ids_exclude{i_field,i_exclude} = ids_EXCLUDE{i_field,i_exclude}(isnan(ids_EXCLUDE{i_field,i_exclude})==0);
    end
    
end


ids_well_suited = ids_WELL_SUITED(isnan(ids_WELL_SUITED)==0);
ids_not_suited = ids_NOT_SUITED(isnan(ids_NOT_SUITED)==0);

RASTER_BASINS = river_points(i_point).refBASIN_depth_dam_crosssection;
dam_x = river_points(i_point).refDAM_dam_coords(1,:)' + (river_points(i_point).refDAM_LLcorner(1)-1) -(river_points(i_point).refBASIN_LLcorner(1)-1);
dam_y = river_points(i_point).refDAM_dam_coords(2,:)' + (river_points(i_point).refDAM_LLcorner(2)-1) -(river_points(i_point).refBASIN_LLcorner(2)-1);
dam_heights = river_points(i_point).dam_heights;

dam_point_area_max = 0;
for i_dam = 1:length(dam_heights_potential)
        RASTER_BASIN = RASTER_BASINS(:,:,i_dam);
        ids_basin{i_dam} = RASTER_ID(isnan(RASTER_BASIN)==0);
        ids_dam{i_dam} = sub2ind(size(RASTER_BASIN),dam_y(dam_heights(i_dam,:)>0),dam_x(dam_heights(i_dam,:)>0));
        ids_basin_dam{i_dam} = [ids_basin{i_dam};ids_dam{i_dam}];
        cells_basin(i_dam) = length(ids_basin_dam{i_dam});
        
        for i_field = 1:length(fields)
            for i_exclude = 1:size(ids_EXCLUDE,2)
                cells_exclude{i_field,i_exclude}(i_dam) = sum(ismember(ids_basin_dam{i_dam},ids_exclude{i_field,i_exclude}));
                
                if cells_exclude{i_field,i_exclude}(i_dam)>0 && cells_exclude{i_field,i_exclude}(i_dam) > info_exclude_basin.(fields{i_field}).exclude_threshold{i_exclude}/(grids.cellsize)^2
                    dam_heights_potential(i_dam) = NaN;
                end
            end
        end
        cells_well_suited(i_dam) = sum(ismember(ids_basin_dam{i_dam},ids_well_suited));
        cells_not_suited(i_dam) = sum(ismember(ids_basin_dam{i_dam},ids_not_suited));
        
        % determine catchment area including all river points in the basin
        if i_dam > dam_point_area_max
            id_river_points = [river_points_all.id_grid];
            id_dam_points = sub2ind(size(grids.dem),...
                [river_points(i_point).refDAM_dam_coords(2,:)+river_points(i_point).refDAM_LLcorner(2)-1],...
                [river_points(i_point).refDAM_dam_coords(1,:)+river_points(i_point).refDAM_LLcorner(1)-1]);
            line_river_points_on_dam = find(ismember(id_river_points,id_dam_points)==1);
            if length(line_river_points_on_dam) > 1
                dam_point_acc_tmp = [];
                for i_river_section = 1:length(river_longitudinal_section)
                    lines_intersect = find(ismember(line_river_points_on_dam,river_longitudinal_section{i_river_section})==1);                    
                    if isempty(lines_intersect)==0
                        dam_point_acc_tmp(end+1) = max([river_points_all(line_river_points_on_dam(lines_intersect)).acc]);
                    end
                    line_river_points_on_dam = setdiff(line_river_points_on_dam,line_river_points_on_dam(lines_intersect));
                    if isempty(line_river_points_on_dam)==1
                        break
                    end
                end
                if length(dam_point_acc_tmp)>1
                    dam_point_area(i_dam,1) = sum(dam_point_acc_tmp)*(grids.cellsize)^2/1e6;
                else
                    dam_point_area(i_dam:length(dam_heights_potential),1) = sum(dam_point_acc_tmp)*(grids.cellsize)^2/1e6;
                    dam_point_area_max = length(dam_heights_potential);
                end
                
            else
                dam_point_area(i_dam:length(dam_heights_potential),1) = river_points_all(i_point_all).acc*(grids.cellsize)^2/1e6;
                dam_point_area_max = length(dam_heights_potential);
            end
        end
        
end

if sum(isnan(dam_heights_potential)==1) == length(dam_heights_potential)
    % too many exclusion cells in the basin
    river_points(i_point).exit_code = 18; 
    river_points_all(i_point_all).exit_code = 18;
    return
end

% specific volume
specific_volume = river_points(i_point).basin_volumes(:,1)./(dam_point_area(:,1).*1e3);
if isempty(sV_min) ==0
    
% specific volume of the basin is smaller than the defined minimum specific volume
[dam_ids_exclude,~] = find(specific_volume < sV_min);
dam_heights_potential(1,dam_ids_exclude) = NaN;
clear dam_ids_exclude

if sum(isnan(dam_heights_potential)==1) == length(dam_heights_potential)
    river_points(i_point).exit_code = 16; % dam is too short
    river_points_all(i_point_all).exit_code = 16; % dam is too short
    return
end

% specific volume of the basin is larger than the defined maximum specific volume
[dam_ids_exclude,~] = find(specific_volume > sV_max);
dam_heights_potential(1,dam_ids_exclude) = NaN;
clear dam_ids_exclude

if sum(isnan(dam_heights_potential)==1) == length(dam_heights_potential)
    river_points(i_point).exit_code = 17; % dam is too short
    river_points_all(i_point_all).exit_code = 17; % dam is too short
    return
end

end

%% calculate evaluation criteria
    % (is required for later comparison among the basins)
    
    dam_ids_exclude = dam_heights_potential;
    dam_ids_exclude(isnan(dam_ids_exclude)==0) = 1;
    
    % criterium 1: dam volume per basin storage volume
    crit1_damvol_per_storage(1,:) = river_points(i_point).dam_volumes(:,1)./river_points(i_point).basin_volumes(:,1);%.*dam_ids_exclude';
    
    % criterium 2: basin area per basin storage volume
    crit2_area_per_storage(1,:) = river_points(i_point).basin_areas(:,1)./river_points(i_point).basin_volumes(:,1);%.*dam_ids_exclude';
    
    % criterium 3: share of well-suited cells
    crit3_share_wellsuited_cells(1,:) = cells_well_suited./cells_basin;
    
    % criterium 4: share of not-suited cells
    crit4_share_notsuited_cells(1,:) = cells_not_suited./cells_basin;
    

%% check if there are more than one dam to evaluate
if length(dam_heights_potential) == 1
    crit1_damvol_per_storage_norm = crit1_damvol_per_storage;
    crit2_area_per_storage_norm = crit2_area_per_storage;
    crit3_share_wellsuited_cells_norm = crit3_share_wellsuited_cells;
    crit4_share_notsuited_cells_norm = crit4_share_notsuited_cells;
    w3_share_well_suited_singleBasin = w3_share_well_suited;
    w4_share_not_suited_singleBasin = w4_share_not_suited;
else      
    %% normalize evaluation criteria
    
    % criterium 1: dam volume per basin storage volume
    crit1_damvol_per_storage_norm = 1-(crit1_damvol_per_storage - min(crit1_damvol_per_storage,[],'omitnan'))/(max(crit1_damvol_per_storage,[],'omitnan') - min(crit1_damvol_per_storage,[],'omitnan'));
    
    % criterium 2: basin area per basin storage volume
    crit2_area_per_storage_norm = 1-(crit2_area_per_storage - min(crit2_area_per_storage,[],'omitnan'))/(max(crit2_area_per_storage,[],'omitnan') - min(crit2_area_per_storage,[],'omitnan'));
    
    % criterium 3: share of well-suited cells
    if max(crit3_share_wellsuited_cells - min(crit3_share_wellsuited_cells,[],'omitnan')) == 0
        w3_share_well_suited_singleBasin = 0;
        crit3_share_wellsuited_cells_norm = zeros(1,size(crit3_share_wellsuited_cells,2));
    else
        w3_share_well_suited_singleBasin = w3_share_well_suited;
        crit3_share_wellsuited_cells_norm = (crit3_share_wellsuited_cells - min(crit3_share_wellsuited_cells,[],'omitnan'))/(max(crit3_share_wellsuited_cells,[],'omitnan') - min(crit3_share_wellsuited_cells,[],'omitnan'));
    end
    
    % criterium 4: share of not-suited cells
    if max(crit4_share_notsuited_cells - min(crit4_share_notsuited_cells,[],'omitnan')) == 0
        w4_share_not_suited_singleBasin = 0;
        crit4_share_notsuited_cells_norm = ones(1,size(crit4_share_notsuited_cells,2));
    else
        w4_share_not_suited_singleBasin = w4_share_not_suited;
        crit4_share_notsuited_cells_norm = 1-(crit4_share_notsuited_cells - min(crit4_share_notsuited_cells,[],'omitnan'))/(max(crit4_share_notsuited_cells,[],'omitnan') - min(crit4_share_notsuited_cells,[],'omitnan'));
    end
end
    
    %% calculation of the objective function and select potential dam
    
    % adapt shares so that there sum will be 1
    w1 = w1_damVolume_per_basinVolume/(w1_damVolume_per_basinVolume+w2_basinArea_per_basinVolume+w3_share_well_suited_singleBasin+w4_share_not_suited_singleBasin);
    w2 = w2_basinArea_per_basinVolume/(w1_damVolume_per_basinVolume+w2_basinArea_per_basinVolume+w3_share_well_suited_singleBasin+w4_share_not_suited_singleBasin);
    w3 = w3_share_well_suited_singleBasin/(w1_damVolume_per_basinVolume+w2_basinArea_per_basinVolume+w3_share_well_suited_singleBasin+w4_share_not_suited_singleBasin);
    w4 = w4_share_not_suited_singleBasin/(w1_damVolume_per_basinVolume+w2_basinArea_per_basinVolume+w3_share_well_suited_singleBasin+w4_share_not_suited_singleBasin);
    
    objective_function = w1 * crit1_damvol_per_storage_norm + ...
        w2 * crit2_area_per_storage_norm + ...
        w3 * crit3_share_wellsuited_cells_norm + ...
        w4 * crit4_share_notsuited_cells_norm;
    
    i_dam = find((objective_function.*dam_ids_exclude)==max(objective_function.*dam_ids_exclude),1,'first');
% end

%% determine global basin ids

[basin_y,basin_x] = ind2sub(size(river_points(i_point).refBASIN_depth_dam_crosssection(:,:,i_dam)),ids_basin{i_dam});
basin_x_global = basin_x + river_points(i_point).refBASIN_LLcorner(1)-1;
basin_y_global = basin_y + river_points(i_point).refBASIN_LLcorner(2)-1;
dam_x_global = dam_x + river_points(i_point).refBASIN_LLcorner(1)-1;
dam_y_global = dam_y + river_points(i_point).refBASIN_LLcorner(2)-1;

ids_basin_global = sub2ind(size(grids.dem),[basin_y_global;dam_y_global],[basin_x_global;dam_x_global]);
ids_basin_global = unique(ids_basin_global);

%% save the selected potential dam

river_points(i_point).refGK_dam_coords(1,:) = ...
    (river_points(i_point).refDAM_dam_coords(1,:)*grids.cellsize - grids.cellsize) + (river_points(i_point).refDAM_LLcorner(1)-1)*grids.cellsize + refGK_LLcornerGLOBAL(1);
river_points(i_point).refGK_dam_coords(2,:) = ...
    refGK_LLcornerGLOBAL(2)+(length(grids.y)-1)*grids.cellsize -(river_points(i_point).refDAM_LLcorner(2)-1)*grids.cellsize-(river_points(i_point).refDAM_dam_coords(2,:)*grids.cellsize-grids.cellsize);
river_points(i_point).refDAM_dam_coords = [];

river_points(i_point).dam_height = dam_heights_potential(isnan(dam_heights_potential)==0)';%(i_dam);
river_points(i_point).dam_top_elev = river_points(i_point).dem + dam_heights_potential(isnan(dam_heights_potential)==0)';
river_points(i_point).dam_row_selected = i_dam - sum(isnan(dam_heights_potential(1:i_dam-1))==1);
river_points(i_point).objective_function_dams = objective_function(isnan(dam_heights_potential)==0)';
river_points(i_point).curve_potential_dam_basin = [river_points(i_point).dam_heights_eval',...
    river_points(i_point).dam_volumes(:,1),river_points(i_point).dam_axis_lengths,...
    river_points(i_point).basin_volumes(:,1),river_points(i_point).basin_areas(:,1),...
    specific_volume,...
    objective_function',dam_ids_exclude',dam_height_max_real,...
    crit1_damvol_per_storage',crit2_area_per_storage',...
    crit3_share_wellsuited_cells',crit4_share_notsuited_cells',...
    cells_not_suited',dam_point_area];
river_points(i_point).dam_volumes = river_points(i_point).dam_volumes(isnan(dam_heights_potential)==0,:);
river_points(i_point).basin_volumes = river_points(i_point).basin_volumes(isnan(dam_heights_potential)==0,:);
river_points(i_point).basin_areas = river_points(i_point).basin_areas(isnan(dam_heights_potential)==0,:);
river_points(i_point).dam_heights = river_points(i_point).dam_heights(isnan(dam_heights_potential)==0,:);
river_points(i_point).dam_elevations = river_points(i_point).dam_elevations(isnan(dam_heights_potential)==0,:);
river_points(i_point).dam_axis_lengths = river_points(i_point).dam_axis_lengths(isnan(dam_heights_potential)==0,:);
river_points(i_point).dam_axis_heights = river_points(i_point).dam_axis_heights(isnan(dam_heights_potential)==0,:);

river_points(i_point).crit1_damvol_per_storage = crit1_damvol_per_storage(isnan(dam_heights_potential)==0)';
river_points(i_point).crit2_area_per_storage = crit2_area_per_storage(isnan(dam_heights_potential)==0)';
river_points(i_point).crit3_share_wellsuited_cells = crit3_share_wellsuited_cells(isnan(dam_heights_potential)==0)';
river_points(i_point).crit4_share_notsuited_cells = crit4_share_notsuited_cells(isnan(dam_heights_potential)==0)';
river_points(i_point).crit_selected = [crit1_damvol_per_storage(1,i_dam),...
    crit2_area_per_storage(1,i_dam),crit3_share_wellsuited_cells(1,i_dam),crit4_share_notsuited_cells(1,i_dam)];

river_points(i_point).ids_basin = ids_basin_global;
river_points(i_point).refBASIN_depth_dam_crosssection = river_points(i_point).refBASIN_depth_dam_crosssection(:,:,i_dam);
river_points(i_point).refBASIN_depth_wall_crosssection = river_points(i_point).refBASIN_depth_wall_crosssection(:,:,i_dam);
