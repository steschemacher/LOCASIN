function fun_8_determine_depth_storage_area_curves(i_point)
% ######  calculate dapth-storage-area-curves  ######
% function to calculate curves of the relation between water depth, storage
% volume and flooding area for the basins of the selected combination
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global discretization_number grids basins_selected

% discretization heights for the impounding levels
heights(:,1) = 0:basins_selected(i_point).dam_height/(discretization_number-1):basins_selected(i_point).dam_height;

% temporary array to avoid loops for the evaluation
heights_temp(1,1,:) = (-1)*(heights-max(heights));

% evaluation for a dam-crosssection and a wall-crosssection
for i_BASIN = 1:2
    if i_BASIN == 1
        BASIN_DEPTH = basins_selected(i_point).refBASIN_depth_dam_crosssection;
    else
        BASIN_DEPTH = basins_selected(i_point).refBASIN_depth_wall_crosssection;
    end
    
    % expand basin depths to the z-direction (one layer per discetization level)
    BASIN_DEPTH = repmat(BASIN_DEPTH,1,1,length(heights));
    
    % exlude depth, which are below the temporary heights (opposed direction)
    BASIN_DEPTH(BASIN_DEPTH<heights_temp)=NaN;
    
    % calculate true basin depth (for the right extent)
    BASIN_DEPTH = BASIN_DEPTH-heights_temp;
    
    % calculate volume and area of the respective impounding levels
    volume(:,1) = squeeze(sum(sum(BASIN_DEPTH,1,'omitnan'),2,'omitnan'))*(grids.cellsize)^2;
    area(:,1) = squeeze(sum(sum(BASIN_DEPTH>=0,1,'omitnan'),2,'omitnan'))*(grids.cellsize)^2;
    
    % combine results
    if i_BASIN == 1
        ShA_dam = [volume,heights,area];
    else
        ShA_wall = [volume,heights,area];
    end
end

% save results to struct
basins_selected(i_point).ShA_dam = ShA_dam;
basins_selected(i_point).ShA_wall = ShA_wall;

if isempty(ShA_dam)==1 || isempty(ShA_wall)==1
    fprintf('Error during the curve calculation.')
end
