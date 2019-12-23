function [dams_selected] = fun_7_selection_of_basin_combination(dam_points)
% ######  select best basin combination  ######
% function to select the best basin combination based on the evaluation
% criteria and overlapping flooding areas

% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global w1_damVolume_per_basinVolume w2_basinArea_per_basinVolume 
global w3_share_well_suited w4_share_not_suited

%% get evaluation criteria

% select criteria
crit_selected = cat(1,dam_points.crit_selected);
crit1_damvol_per_storage = crit_selected(:,1);
crit2_area_per_storage = crit_selected(:,2);
crit3_share_wellsuited_cells = crit_selected(:,3);
crit4_share_notsuited_cells = crit_selected(:,4);


%% normalize evaluation criteria

% criterium 1: dam volume per basin storage volume
crit1_damvol_per_storage_norm = 1-(crit1_damvol_per_storage - min(crit1_damvol_per_storage,[],'omitnan'))/max(crit1_damvol_per_storage - min(crit1_damvol_per_storage,[],'omitnan'));

% criterium 2: basin area per basin storage volume
crit2_area_per_storage_norm = 1-(crit2_area_per_storage - min(crit2_area_per_storage,[],'omitnan'))/max(crit2_area_per_storage - min(crit2_area_per_storage,[],'omitnan'));

% criterium 3: share of well-suited cells
if max(crit3_share_wellsuited_cells - min(crit3_share_wellsuited_cells,[],'omitnan')) == 0
    w3_share_well_suited = 0;
    crit3_share_wellsuited_cells_norm = zeros(1,size(crit3_share_wellsuited_cells,2));
else
    crit3_share_wellsuited_cells_norm = (crit3_share_wellsuited_cells - min(crit3_share_wellsuited_cells,[],'omitnan'))/max(crit3_share_wellsuited_cells - min(crit3_share_wellsuited_cells,[],'omitnan'));
end

% criterium 4: share of not-suited cells
if max(crit4_share_notsuited_cells - min(crit4_share_notsuited_cells,[],'omitnan')) == 0
    w4_share_not_suited = 0;
    crit4_share_notsuited_cells_norm = ones(1,size(crit4_share_notsuited_cells,2));
else
    crit4_share_notsuited_cells_norm = 1-(crit4_share_notsuited_cells - min(crit4_share_notsuited_cells,[],'omitnan'))/max(crit4_share_notsuited_cells - min(crit4_share_notsuited_cells,[],'omitnan'));
end


%% calculation of the objective function and select potential dam

% adapt shares so that there sum will be 1
w1 = w1_damVolume_per_basinVolume/(w1_damVolume_per_basinVolume+w2_basinArea_per_basinVolume+w3_share_well_suited+w4_share_not_suited);
w2 = w2_basinArea_per_basinVolume/(w1_damVolume_per_basinVolume+w2_basinArea_per_basinVolume+w3_share_well_suited+w4_share_not_suited);
w3 = w3_share_well_suited/(w1_damVolume_per_basinVolume+w2_basinArea_per_basinVolume+w3_share_well_suited+w4_share_not_suited);
w4 = w4_share_not_suited/(w1_damVolume_per_basinVolume+w2_basinArea_per_basinVolume+w3_share_well_suited+w4_share_not_suited);

% calculation of the objective function
objective_function_cell = num2cell(w1 * crit1_damvol_per_storage_norm + ...
    w2 * crit2_area_per_storage_norm + ...
    w3 * crit3_share_wellsuited_cells_norm + ...
    w4 * crit4_share_notsuited_cells_norm);

[dam_points.objective_function] = objective_function_cell{:};

%% Select the best basins

%%% sort basins by the objective function:
% good values (large) on top
% bad values (small) at the bottom
[~,index] = sortrows([dam_points.objective_function].',-1);
dam_points = dam_points(index);
clear index

becken = [dam_points.id];
dam_points_tmp = dam_points;

% choose the best basin, if basin points are overlapping
% (best is on the top, which means it is the first to choose)
m=1;
while isempty(becken)==0
    
    % include all potential basins
    becken_ids = [dam_points([dam_points.id]==becken(1)).ids_basin];
    
    % choose the first basin for the combination (which is the best)
    dams_selected(m) = dam_points([dam_points.id] == becken(1));
    m = m+1;
    
    % delete all basins which are overlapping with the chosen one from
    % the list
    for i = 1:length(dam_points_tmp)
        if sum(ismember([dam_points_tmp(i).ids_basin],becken_ids))>=1
            % exclude overlapping basins from the basin list
            becken = becken(becken~=dam_points_tmp(i).id);
        end
    end
    clear dam_points_tmp
    % update the dam points which are left for the analysis
    dam_points_tmp = dam_points(ismember([dam_points.id],becken));
    
end
clear dam_points_tmp

