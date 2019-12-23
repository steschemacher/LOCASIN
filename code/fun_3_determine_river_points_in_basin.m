function fun_3_determine_river_points_in_basin(i_point,i_point_all)
% ######  estimate basin extent  ######
% function to determine all potential points in the basin to extimate the
% required extent of the raster for the basin calculation
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019


global river_points river_points_all debug_on dam_height_min
global river_longitudinal_section grids point_boundaries

basin_points_elevation_max = river_points(i_point).dam_top_elev;
basin_points = [];
endpoints = [];

% the dam has the favored height starting condition
river_points(i_point).backwater_outside_dem = 0;

for i_start = 1:length(river_longitudinal_section)
    endpunkt = [];
    % check if river point is in the river section
    loc_river_point = find(river_longitudinal_section{i_start}==river_points(i_point).id);
    if isempty(loc_river_point)==0 && isnan(loc_river_point)==0
        % check if the end of the river section is within the global DEM
        potential_river_points_DEM = [river_points_all(river_longitudinal_section{i_start}(1:loc_river_point)).dem];
        loc_point_above_crest = find(potential_river_points_DEM...
            >=river_points(i_point).dam_top_elev,1,'last');
        
        if isempty(loc_point_above_crest)==1 || isnan(loc_point_above_crest)==1      
            loc_point_above_crest = 1;
            id_point_river_all = find([river_points_all.id]==river_longitudinal_section{i_start}(1));
            if river_points_all(id_point_river_all).col<=2 || river_points_all(id_point_river_all).col>=size(grids.dem,2)-1 ||...
                    river_points_all(id_point_river_all).row<=2 || river_points_all(id_point_river_all).row>=size(grids.dem,1)-1
                % determine the highest upstream neighbor
                basin_points_elevation_max = min(max(potential_river_points_DEM),basin_points_elevation_max);            
            end
        end
        
        % add points to basin_points
        basin_points = [basin_points,...
            fliplr(river_longitudinal_section{i_start}(loc_point_above_crest+1:loc_river_point-1))];
        endpunkt = river_longitudinal_section{i_start}(loc_point_above_crest);
        
        % check if dam endpoints exist for the uppermost neighbor in the
        % basin (or if it had been excluded from the analysis due to the
        % computational effort)
        if river_points_all(endpunkt).exit_code == 5
            p_river_point = find([river_points.id]==endpunkt);
            p_river_point_all = find([river_points_all.id]==endpunkt);
            if isempty(p_river_point)==0
                fun_2_determine_shortest_dam(p_river_point,p_river_point_all)
            end
        end
        
    end
    endpoints = [endpoints,endpunkt];
end

dam_height = floor((basin_points_elevation_max - river_points(i_point).dem)*10)/10;

if dam_height >= dam_height_min
    basin_points = unique(basin_points,'stable');
    ids_inBecken = ismember([river_points_all.id],basin_points);
    basin_points_endpoint_coords = cat(2,river_points_all(ids_inBecken==1).refGLOBAL_damEndpoints_coords);
    
    % determine points within the basin
    % (river points upstream + dam endpoints of these river points)
    % ids_inBecken = ismember([river_points.id],[river_points(i_point).basin_points]);
    basin_basinpoints = [river_points_all(i_point_all).refGLOBAL_damEndpoints_coords,...
        basin_points_endpoint_coords,...
        [[river_points_all(basin_points).col];[river_points_all(basin_points).row]]];
    
    
    % select part of the DEM --> refBASIN DEM
    Ymin = max(1,round(min(basin_basinpoints(2,:))-point_boundaries/4));
    Ymax = min(size(grids.dem,1),round(max(basin_basinpoints(2,:))+point_boundaries/4));
    Xmin = max(1,round(min(basin_basinpoints(1,:))-point_boundaries/4));
    Xmax = min(size(grids.dem,2),round(max(basin_basinpoints(1,:))+point_boundaries/4));
    
    if dam_height < river_points(i_point).dam_height
        % basin end point exists (river point), but dam height had to be
        % reduced due to a missing basin end point for the original dam height
        river_points(i_point).exit_code = -2;
        river_points_all(i_point_all).exit_code = -2;
    end
    % save results
    river_points(i_point).dam_height = dam_height;
    river_points(i_point).dam_top_elev = river_points(i_point).dem + dam_height;
    river_points(i_point).refGLOBAL_basin_boundaries = [Xmin,Xmax;Ymin,Ymax];
    river_points(i_point).basin_points = basin_points;   % river points in the basin (ids)
    river_points(i_point).endpunkt = endpoints;          % most upstream river point(s) in the basin (ids)
    
else    
    river_points(i_point).exit_code = 7;  % no basin end point exists (river point)
    river_points_all(i_point_all).exit_code = 7;  % no basin end point exists (river point)    
end


if debug_on == 1
    beckenpunkt = basin_points;
	% definition of colormap
    anzahl = 10;
    orange = [177, 108, 37]/255;
    rot = [141, 36, 37]/255;
    blau1 = [105, 151, 201]/255;
    blau4 = [36, 68, 101]/255;
    orange_rot = [linspace(orange(1),rot(1),anzahl)',linspace(orange(2),rot(2),anzahl)',linspace(orange(3),rot(3),anzahl)'];
    blau14 = [linspace(blau1(1),blau4(1),anzahl)',linspace(blau1(2),blau4(2),anzahl)',linspace(blau1(3),blau4(3),anzahl)'];
    cmap = [flipud(blau14);orange_rot];
    
    % coordinates of river points in the basin
    river_points_in_basin = [[river_points_all(beckenpunkt).col];[river_points_all(beckenpunkt).row]];
        
    % selection of the grid (global water depth)
    GRID =  grids.dem - (river_points(i_point).dem + river_points(i_point).dam_height);
    
    % definition of the range of the colormap (based on local water depth grid)
    max_abs_val = river_points(i_point).dam_height;
    
    % plot data
    figure('color',[1,1,1],'position',[1,1,1000,1000])
    hold on
    % plot grid with colorbar
    imagesc(GRID)
    caxis([-max_abs_val,max_abs_val])
    colormap(cmap)
    cb = colorbar;
    ylabel(cb,'elevation relative to the elevation of the dam top (= water depth)')
    % plot points of dam and river points in the basin
    p0 = plot(river_points_in_basin(1,:),river_points_in_basin(2,:),'.k');
    p1 = plot(river_points(i_point).col,river_points(i_point).row,'*r');
    p2 = plot(river_points(i_point).refGLOBAL_damEndpoints_coords(1,:),river_points(i_point).refGLOBAL_damEndpoints_coords(2,:),'or');
    p3 = plot(river_points(i_point).refGLOBAL_damBoundarypoints_coords(1,:),river_points(i_point).refGLOBAL_damBoundarypoints_coords(2,:),'--og','markersize',8);
    % legend
    legend([p0,p1,p2,p3],'river points in the basin','river point of the dam','dam endpoints','points at the refDAM-boundaries')
    axis equal tight
    set(gca,'ydir','reverse')
    
    keyboard
    
end
