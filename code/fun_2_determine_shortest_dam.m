function fun_2_determine_shortest_dam(i_point,i_point_all)
% ######  shortest dam orientation  ######
% function to determine the dam orientation with the shortest dam length
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019


global debug_on point_boundaries river_points river_points_all 
global grids dam_length_max
global dam_dist_eval dam_height_max dam_height_min
global neighbors_exclude_distance

n_tries_max = 10;
shift_border = round(point_boundaries/(n_tries_max/5));

% select relevant part of the DEM for the dam
% (cut dem around potential dam point for the analysis)
Xmin = ceil(river_points(i_point).col-dam_length_max/grids.cellsize/3);
Xmax = floor(river_points(i_point).col+dam_length_max/grids.cellsize/3);
Ymin = ceil(river_points(i_point).row-dam_length_max/grids.cellsize/3);
Ymax = floor(river_points(i_point).row+dam_length_max/grids.cellsize/3);
refDAM_dem = grids.dem(Ymin:Ymax,Xmin:Xmax);

if isnan(sum(refDAM_dem,'all'))==1
    % there are NaN values in the analysis plot
    % --> the river point is too close to the data border
    river_points(i_point).exit_code = 1;
    river_points_all(i_point_all).exit_code = 1;
    return
end

point_xy = [river_points(i_point).col;river_points(i_point).row]-[Xmin;Ymin]+[1;1];


% determine potential inundation area 
% (dam endpoint needs to be above water level)
refDAM_above_waterlevel = refDAM_dem - (river_points(i_point).dem + dam_height_max);
refDAM_above_waterlevel(refDAM_above_waterlevel>=0)=1;
refDAM_above_waterlevel(refDAM_above_waterlevel<0)=0;

% determine extent of the inundation area
Xmin_water = find(sum(refDAM_above_waterlevel,1,'omitnan')<size(refDAM_above_waterlevel,1),1,'first');
Xmax_water = find(sum(refDAM_above_waterlevel,1,'omitnan')<size(refDAM_above_waterlevel,1),1,'last');
Ymin_water = find(sum(refDAM_above_waterlevel,2,'omitnan')<size(refDAM_above_waterlevel,2),1,'first');
Ymax_water = find(sum(refDAM_above_waterlevel,2,'omitnan')<size(refDAM_above_waterlevel,2),1,'last');

% fit the clipped DEM in order to reduce land area
Xmax = Xmin + min(size(refDAM_above_waterlevel,1),Xmax_water+2)-1;
Xmin = Xmin + max(1,Xmin_water-2)-1;
Ymax = Ymin + min(size(refDAM_above_waterlevel,2),Ymax_water+2)-1;
Ymin = Ymin + max(1,Ymin_water-2)-1;

% reclip the DEM
refDAM_dem = grids.dem(Ymin:Ymax,Xmin:Xmax);
point_xy = [river_points(i_point).col;river_points(i_point).row]-[Xmin;Ymin]+[1;1];

no_tries_necessary = zeros(1,8);
for n_tries = 1:n_tries_max
    
    if sum(no_tries_necessary)==8
        break
    end  
    
    % determine the new inundation area
    refDAM_above_waterlevel = refDAM_dem - (river_points(i_point).dem + dam_height_max);
    refDAM_above_waterlevel(refDAM_above_waterlevel>=0)=1;
    refDAM_above_waterlevel(refDAM_above_waterlevel<0)=0;
    
    % enlarge boundaries if necessary
    % (Boundaries need to be enlarged, if no valley exits in the clipped
    % area. This is usually indicated, if there is only water at one boundary line
    % (== 0) and at the same time not one "only-water"-line in the
    % perpendicular direction. (there could be a narrow valley downstream
    % or upstream of the river point, which would also give that result,
    % but this region would be much better than the investigated. Thus, the
    % investigated river point can be excluded)
    % left border
    if (sum(refDAM_above_waterlevel(:,1),'omitnan')==0 && min(sum(refDAM_above_waterlevel,2,'omitnan'))~=0) ...
            || sum(refDAM_above_waterlevel,'all')==0
        Xmin = max(1,Xmin-shift_border);
    else
        no_tries_necessary(1)=1;
    end
    % right border
    if (sum(refDAM_above_waterlevel(:,end),'omitnan')==0 && min(sum(refDAM_above_waterlevel,2,'omitnan'))~=0) ...
            || sum(refDAM_above_waterlevel,'all')==0
%         Xmin = max(1,Xmin-shift_border);
        Xmax = min(length(grids.x),Xmax+shift_border);
    else
        no_tries_necessary(2)=1;
    end
    % bottom border (Ymax)
    if (sum(refDAM_above_waterlevel(1,:),'omitnan')==0 && min(sum(refDAM_above_waterlevel,1,'omitnan'))~=0) ...
            || sum(refDAM_above_waterlevel,'all')==0
%         Xmin = max(1,Xmin-shift_border);
        Ymin = max(1,Ymin-shift_border);
    else
        no_tries_necessary(3)=1;
    end
    % top border (Ymin)
    if (sum(refDAM_above_waterlevel(end,:),'omitnan')==0 && min(sum(refDAM_above_waterlevel,1,'omitnan'))~=0) ...
            || sum(refDAM_above_waterlevel,'all')==0
%         Xmin = max(1,Xmin-shift_border);
        Ymax = min(length(grids.y),Ymax+shift_border);
    else
        no_tries_necessary(4)=1;
    end
    % top-left border
    if sum(refDAM_above_waterlevel(end,:),'omitnan')==0 && sum(refDAM_above_waterlevel(:,1),'omitnan')==0
        Xmin = max(1,Xmin-shift_border);
        Ymax = min(length(grids.y),Ymax+shift_border);
    else
        no_tries_necessary(5)=1;
    end
    % top-right border
    if sum(refDAM_above_waterlevel(end,:),'omitnan')==0 && sum(refDAM_above_waterlevel(:,end),'omitnan')==0
        Xmax = min(length(grids.x),Xmax+shift_border);
        Ymax = min(length(grids.y),Ymax+shift_border);
    else
        no_tries_necessary(6)=1;
    end
    % bottom-right border
    if sum(refDAM_above_waterlevel(1,:),'omitnan')==0 && sum(refDAM_above_waterlevel(:,end),'omitnan')==0
        Xmax = min(length(grids.x),Xmax+shift_border);
        Ymin = max(1,Ymin-shift_border);
    else
        no_tries_necessary(7)=1;
    end
    % bottom-left border
    if sum(refDAM_above_waterlevel(1,:),'omitnan')==0 && sum(refDAM_above_waterlevel(:,1),'omitnan')==0
        Xmin = max(1,Xmin-shift_border);
        Ymin = max(1,Ymin-shift_border);
    else
        no_tries_necessary(8)=1;
    end
    refDAM_dem = grids.dem(Ymin:Ymax,Xmin:Xmax);
end

point_xy = [river_points(i_point).col;river_points(i_point).row]-[Xmin;Ymin]+[1;1];
refDAM_size = [size(refDAM_dem,2),size(refDAM_dem,1)];


for dam_height = dam_height_max:-dam_dist_eval:dam_height_min %sort(dam_height_max,'descend')
    % determine the maximum water depth for the selected dem around the point
    refDAM_waterdepth_max = refDAM_dem - (river_points(i_point).dem + dam_height);
    
    % dam endpoint needs to be above water level
    refDAM_above_waterlevel = refDAM_waterdepth_max;
    refDAM_above_waterlevel(refDAM_above_waterlevel>=0)=1;
    refDAM_above_waterlevel(refDAM_above_waterlevel<0)=0;
    point_xy = [river_points(i_point).col;river_points(i_point).row]-[Xmin;Ymin]+[1;1];
    
    % the DEM has to be a valley 
    if sum(no_tries_necessary)<8
        % --> 3 of the edges have to be land (= 1) or 2 edges on opposite sides
        if sum([refDAM_above_waterlevel(1,1),refDAM_above_waterlevel(1,end),...
                refDAM_above_waterlevel(end,1),refDAM_above_waterlevel(end,end)],'omitnan') < 3
            if sum([refDAM_above_waterlevel(1,1),refDAM_above_waterlevel(end,end)],'omitnan')<2 ...
                    && sum([refDAM_above_waterlevel(end,1),refDAM_above_waterlevel(1,end)],'omitnan')<2
                if dam_height > dam_height_max
                    continue
                else
                    river_points(i_point).exit_code = 6;  % no dam of this height was found for the point
                    river_points_all(i_point_all).exit_code = 6;  % no dam of this height was found for the point
                    return
                end
            end
        end
        % --> 2 border sides have to have water (= 0)
        sum_edges = zeros(1,4);
        if sum(refDAM_above_waterlevel(1,:),'omitnan')<size(refDAM_above_waterlevel,2)
            sum_edges(1) = 1;
        end
        if sum(refDAM_above_waterlevel(end,:),'omitnan')<size(refDAM_above_waterlevel,2)
            sum_edges(2) = 1;
        end
        if sum(refDAM_above_waterlevel(:,1),'omitnan')<size(refDAM_above_waterlevel,1)
            sum_edges(3) = 1;
        end
        if sum(refDAM_above_waterlevel(:,end),'omitnan')<size(refDAM_above_waterlevel,1)
            sum_edges(4) = 1;
        end
        if sum(sum_edges) < 2
            if dam_height > dam_height_max
                continue
            else
                river_points(i_point).exit_code = 6;  % no dam of this height was found for the point
                river_points_all(i_point_all).exit_code = 6;  % no dam of this height was found for the point
                return
            end
        end
    end
    refDAM_above_waterlevel(refDAM_above_waterlevel==0)=NaN;
    
    
    % determine the points at which the dam ends at the valley borders
    % (minimum dam width)
    clear points min_dist_dam mat_num
    
    %% define thresholds
    
    winkel_lower_bound = pi/2 - .01;%1.517;%
    winkel_upper_bound = pi/2 + .01;%1.597;%
    
    %% definition of rasters for the dam analysis
    % distance in x- and y-direction to the river point
    % (the river point is in the middle of the clipped DEM)
    refDAM_distx = repmat(-point_xy(1)+1:refDAM_size(1)-point_xy(1),refDAM_size(2),1);
    refDAM_disty = repmat([-point_xy(2)+1:refDAM_size(2)-point_xy(2)]',1,refDAM_size(1));
    
    % distance of the cells from the center (= river point)
    refDAM_dist = sqrt(refDAM_distx.^2+refDAM_disty.^2);
    
    if debug_on == 1
        figure('color',[1,1,1],'position',[1,1,1500,500])
        subplot(1,3,1)
        img = imagesc(refDAM_distx);
        set(img,'AlphaData',~isnan(refDAM_distx))
        colorbar
        axis equal tight
        title('refDAM_distx','interpreter','none')
        subplot(1,3,2)
        img = imagesc(refDAM_disty);
        set(img,'AlphaData',~isnan(refDAM_disty))
        colorbar
        axis equal tight
        title('refDAM_disty','interpreter','none')
        subplot(1,3,3)
        img = imagesc(refDAM_dist);
        set(img,'AlphaData',~isnan(refDAM_dist))
        colorbar
        axis equal tight
        title('refDAM_dist','interpreter','none')
    end
    
    % angles to describe the directions of the dam segments (left and right of river point)
    refDAM_anglA = atan(refDAM_disty./refDAM_distx); % Angle A(alpha) is between horizontal line and direction-to-point-line
    refDAM_anglB = atan(refDAM_distx./refDAM_disty); % Angle B(beta) is between vertical line and direction-to-point-line
    
    
    % raster to describe if dam endpoint is in the left part (=1) or the right part (=2)
    raster_location12 = zeros(numel(refDAM_dem),1);
    index_dampoint = sub2ind(size(refDAM_dem),point_xy(2),point_xy(1));
    raster_location12(1:index_dampoint)=1;
    raster_location12(index_dampoint+1:end)=2;
    
    refDAM_loc12 = reshape(raster_location12,size(refDAM_anglA));
    
    if debug_on == 1
        figure('color',[1,1,1],'position',[1,1,1500,500])
        subplot(1,3,1)
        img = imagesc(refDAM_anglA);
        set(img,'AlphaData',~isnan(refDAM_anglA))
        colorbar
        axis equal tight
        title('refDAM_anglA','interpreter','none')
        subplot(1,3,2)
        img = imagesc(refDAM_anglB);
        set(img,'AlphaData',~isnan(refDAM_anglB))
        colorbar
        axis equal tight
        title('refDAM_anglB','interpreter','none')
        subplot(1,3,3)
        img = imagesc(refDAM_loc12);
        set(img,'AlphaData',~isnan(refDAM_loc12))
        colorbar
        axis equal tight
        title('refDAM_loc12','interpreter','none')
    end
    
    
    %% reduction of the number of possible dam points
    
    %%% dam endpoint needs to be at the water line
    % (multiply with neighbors to enlarge the water area by one pixel, then
    % substract the result from the original "refDAM_above_waterlevel" to gain
    % the waterline grid)
    refDAM_above_waterlevel_shift = refDAM_above_waterlevel;
    refDAM_above_waterlevel_shift = refDAM_above_waterlevel.*(refDAM_above_waterlevel([2:end,end],:));
    refDAM_above_waterlevel_shift(2:end,:) = refDAM_above_waterlevel_shift(2:end,:).*(refDAM_above_waterlevel([1:end-1],:));
    refDAM_above_waterlevel_shift = refDAM_above_waterlevel_shift.*(refDAM_above_waterlevel(:,[2:end,end]));
    refDAM_above_waterlevel_shift(:,2:end) = refDAM_above_waterlevel_shift(:,2:end).*(refDAM_above_waterlevel(:,[1:end-1]));
    refDAM_waterline = refDAM_above_waterlevel;
    refDAM_waterline(isnan(refDAM_above_waterlevel_shift)==0)=NaN;
    
    
    %%% delete points of the waterline, which have no dry cells behind at the
    % grid boundaries (e.g. islands, where water would flow around
    % --> determine coordinates of the borders behind all raster cells
    % --> select border-coordinates for the shoreline point
    % --> check if there is "land" at these coordinates
    
    % define diagonals from the edges to the river point (every diagonal
    % quarter has either the same maximum x or y coordinate)
    [Y(:,1),X(:,1)] = ind2sub(size(refDAM_above_waterlevel),1:numel(refDAM_above_waterlevel));
    
    
    y1u = (X.*(point_xy(2)-1)/(point_xy(1)-1))...
        + (point_xy(2)-point_xy(1)*((point_xy(2)-1)/(point_xy(1)-1)));
    y1o = (X.*(size(refDAM_above_waterlevel,1)-point_xy(2))/(1-point_xy(1))) ...
        + (point_xy(2)-point_xy(1)*((size(refDAM_above_waterlevel,1)-point_xy(2))/(1-point_xy(1))));
    y2u = (X.*(point_xy(2)-1)/(point_xy(1)-size(refDAM_above_waterlevel,2)) ...
        + (point_xy(2)-point_xy(1)*(point_xy(2)-1)/(point_xy(1)-size(refDAM_above_waterlevel,2))));
    y2o = (X.*(size(refDAM_above_waterlevel,1)-point_xy(2))/(size(refDAM_above_waterlevel,2)-point_xy(1)) ...
        + (point_xy(2)-point_xy(1)*((size(refDAM_above_waterlevel,1)-point_xy(2))/(size(refDAM_above_waterlevel,2)-point_xy(1)))));
    
    % define empty boundary coordinates
    X2 = zeros(size(Y,1),1);
    Y2 = zeros(size(Y,1),1);
    
    % X-coordinates in the left triangular = 1 (betreen y1o and y1u)
    X2(Y>=y1u & Y<=y1o) = 1;    
    % X-coordinates in the right triangular = size(grid,2) (between y2o und y2u)
    X2(Y>=y2u & Y<=y2o) = size(refDAM_above_waterlevel,2);
    % Y-coordinates in the upper triangular = 1 (below y1u and y2u)
    Y2(Y<=y1u & Y<=y2u) = 1; 
    % Y-coordinates in the upper triangular = size(grid,1) (above y1o and y2o)
    Y2(Y>=y1o & Y>=y2o) = size(refDAM_above_waterlevel,1);
    
     % find missing coordinates (X and Y), which are not at the borders
    % (Strahlensatz / intercept theorem)
    mitte_x = repmat(point_xy(1),length(X),1); 
    mitte_y = repmat(point_xy(2),length(X),1);
    index_x = find((Y<y1u & Y<y2u) | (Y>y1o & Y>y2o));
    index_y = find((Y>y1u & Y<y1o) | (Y>y2u & Y<y2o));
    X2(index_x,1) = (mitte_x(index_x)- ((mitte_y(index_x)-Y2(index_x)).*(mitte_x(index_x)-X(index_x,1)))...
        ./(mitte_y(index_x)-Y(index_x,1))); 
    Y2(index_y,1) = (mitte_y(index_y)- ((mitte_x(index_y)-X2(index_y)).*(mitte_y(index_y)-Y(index_y,1)))...
        ./(mitte_x(index_y)-X(index_y,1)));
    X2 = round(X2);
    Y2 = round(Y2);
    
    
    % if the shoreline is on the river point (theoretical consideration for
    % plotting reasons --> will not occur in reality)   
    X2(index_dampoint)=1;
    Y2(index_dampoint)=1;

    
    % generate matrix data
    matX2 = reshape(X2,size(refDAM_above_waterlevel,1),size(refDAM_above_waterlevel,2));
    matY2 = reshape(Y2,size(refDAM_above_waterlevel,1),size(refDAM_above_waterlevel,2));
    
    
    index = sub2ind(size(refDAM_above_waterlevel),Y2,X2);
    index(isnan(index)==1) = 1;
    
    
    % raster with ones, when there is a dry cell at the boundary in this
    % direction
    possible_boundary = refDAM_above_waterlevel(index);
    refDAM_POSSIBLE_BOUNDARY = reshape(possible_boundary,size(refDAM_above_waterlevel));
    
    % reduce the waterline points to the points where dry cells are at the
    % boundary
    refDAM_waterline_boundary = refDAM_waterline .* refDAM_POSSIBLE_BOUNDARY;
    
    
    if debug_on == 1
        figure('color',[1,1,1],'position',[1,1,2000,500])
        subplot(1,4,1)
        img = imagesc(refDAM_above_waterlevel);
        set(img,'AlphaData',~isnan(refDAM_above_waterlevel))
        colorbar
        axis equal tight
        title('refDAM_above_waterlevel','interpreter','none')
        subplot(1,4,2)
        img = imagesc(refDAM_waterline);
        set(img,'AlphaData',~isnan(refDAM_waterline))
        colorbar
        axis equal tight
        title('refDAM_waterline','interpreter','none')
        subplot(1,4,3)
        img = imagesc(refDAM_POSSIBLE_BOUNDARY);
        set(img,'AlphaData',~isnan(refDAM_POSSIBLE_BOUNDARY))
        colorbar
        axis equal tight
        title('refDAM_POSSIBLE_BOUNDARY','interpreter','none')
        subplot(1,4,4)
        img = imagesc(refDAM_waterline_boundary);
        set(img,'AlphaData',~isnan(refDAM_waterline_boundary))
        colorbar
        axis equal tight
        title('refDAM_waterline_boundary','interpreter','none')
    end
    
    
    %% combine the dam analysis rasters (part 1) and the waterline (part 2)
    
    % reduce location12 and distances
    refDAM_POSSIBLE_loc12 = refDAM_loc12.*refDAM_above_waterlevel;
    refDAM_POSSIBLE_DIST =refDAM_dist.*refDAM_above_waterlevel;
    
    % reduce angles alpha and beta
    refDAM_POSSIBLE_ALPHA = refDAM_anglA.*refDAM_waterline_boundary;
    refDAM_POSSIBLE_BETA = refDAM_anglB.*refDAM_waterline_boundary;
    
    % exclude exact horizontal and vertical lines (computational reasons
    refDAM_POSSIBLE_ALPHA(point_xy(2),:) = NaN;
    refDAM_POSSIBLE_ALPHA(:,point_xy(1)) = NaN;
    refDAM_POSSIBLE_BETA(point_xy(2),:) = NaN;
    refDAM_POSSIBLE_BETA(:,point_xy(1)) = NaN;
    refDAM_POSSIBLE_loc12(point_xy(2),:) = NaN;
    refDAM_POSSIBLE_loc12(:,point_xy(1)) = NaN;
    refDAM_waterline_boundary(point_xy(2),:) = NaN;
    refDAM_waterline_boundary(:,point_xy(1)) = NaN;
    refDAM_POSSIBLE_DIST(point_xy(2),:) = NaN;
    refDAM_POSSIBLE_DIST(:,point_xy(1)) = NaN;
    
    if debug_on == 1
        figure('color',[1,1,1],'position',[1,1,2000,500])
        subplot(1,5,1)
        img = imagesc(refDAM_POSSIBLE_ALPHA);
        set(img,'AlphaData',~isnan(refDAM_POSSIBLE_ALPHA))
        colorbar
        axis equal tight
        title('refDAM_POSSIBLE_ALPHA','interpreter','none')
        subplot(1,5,2)
        img = imagesc(refDAM_POSSIBLE_BETA);
        set(img,'AlphaData',~isnan(refDAM_POSSIBLE_BETA))
        colorbar
        axis equal tight
        title('refDAM_POSSIBLE_BETA','interpreter','none')
        subplot(1,5,3)
        img = imagesc(refDAM_POSSIBLE_loc12);
        set(img,'AlphaData',~isnan(refDAM_POSSIBLE_loc12))
        colorbar
        axis equal tight
        title('refDAM_POSSIBLE_loc12','interpreter','none')
        subplot(1,5,4)
        img = imagesc(refDAM_waterline_boundary);
        set(img,'AlphaData',~isnan(refDAM_waterline_boundary))
        colorbar
        axis equal tight
        title('refDAM_waterline_boundary','interpreter','none')
        subplot(1,5,5)
        img = imagesc(refDAM_POSSIBLE_DIST);
        set(img,'AlphaData',~isnan(refDAM_POSSIBLE_DIST))
        colorbar
        axis equal tight
        title('refDAM_POSSIBLE_DIST','interpreter','none')
    end
    
    % save raster data in to array
    array_raster_data(:,1) = 1:numel(refDAM_above_waterlevel);
    array_raster_data(:,2) = reshape(refDAM_POSSIBLE_ALPHA,numel(refDAM_above_waterlevel),1);
    array_raster_data(:,3) = reshape(refDAM_POSSIBLE_BETA,numel(refDAM_above_waterlevel),1);
    array_raster_data(:,4) = reshape(refDAM_POSSIBLE_loc12,numel(refDAM_above_waterlevel),1);
    array_raster_data(:,5) = reshape(refDAM_waterline_boundary,numel(refDAM_above_waterlevel),1);
    array_raster_data(:,6) = reshape(refDAM_POSSIBLE_DIST,numel(refDAM_above_waterlevel),1);
    
    
    %% combine data to determine shortest dam
    
    % reduce array to possible dam
    array_raster_data_small = array_raster_data(array_raster_data(:,5)==1,:);
    
    % separate array in left and right part of the raster (loc12=1, loc12=2)
    array_raster_data_loc1 = array_raster_data_small(array_raster_data_small(:,4)==1,:);
    array_raster_data_loc2 = array_raster_data_small(array_raster_data_small(:,4)==2,:);
    
    
    %%% combine data in matrices for loc1 and loc2
    
    % save angles of alpha in loc1 and beta in loc2 as matrix
    mat_angleA_loc1 = repmat(array_raster_data_loc1(:,2)',[size(array_raster_data_loc2,1),1]);
    mat_angleB_loc2 = repmat(array_raster_data_loc2(:,3),[1,size(array_raster_data_loc1,1)]);
    
    % add angles (preparation to determine, if dam is straight
    mat_angle_sum = (mat_angleA_loc1) + (mat_angleB_loc2);
    
    % save distances in matrix for the combination
    mat_dist_loc1 = repmat(array_raster_data_loc1(:,6)',[size(array_raster_data_loc2,1),1]);
    mat_dist_loc2 = repmat(array_raster_data_loc2(:,6),[1,size(array_raster_data_loc1,1)]);
    
    % determine of the distances of all point-combinations from loc1 and loc2
    mat_dist_dam_all = mat_dist_loc1 + mat_dist_loc2;
    
    if debug_on == 1
        figure('color',[1,1,1],'position',[1,1,2000,1000])
        subplot(2,3,1)
        img = imagesc(mat_angleA_loc1);
        set(img,'AlphaData',~isnan(mat_angleA_loc1))
        colorbar
        axis equal tight
        title('mat_angleA_loc1','interpreter','none')
        subplot(2,3,2)
        img = imagesc(mat_angleB_loc2);
        set(img,'AlphaData',~isnan(mat_angleB_loc2))
        colorbar
        axis equal tight
        title('mat_angleB_loc2','interpreter','none')
        subplot(2,3,3)
        img = imagesc(mat_angle_sum);
        set(img,'AlphaData',~isnan(mat_angle_sum))
        colorbar
        axis equal tight
        title('mat_angle_sum','interpreter','none')
        subplot(2,3,4)
        img = imagesc(mat_dist_loc1);
        set(img,'AlphaData',~isnan(mat_dist_loc1))
        colorbar
        axis equal tight
        title('mat_dist_loc1','interpreter','none')
        subplot(2,3,5)
        img = imagesc(mat_dist_loc2);
        set(img,'AlphaData',~isnan(mat_dist_loc2))
        colorbar
        axis equal tight
        title('mat_dist_loc2','interpreter','none')
        subplot(2,3,6)
        img = imagesc(mat_dist_dam_all);
        set(img,'AlphaData',~isnan(mat_dist_dam_all))
        colorbar
        axis equal tight
        title('mat_dist_dam_all','interpreter','none')
    end
    
    
    %%% combine the angles and select possible combinations
    
    % multiply angles to check if they are in 90° or 180° to each other
    mat_angle_product = mat_angleA_loc1.*mat_angleB_loc2;
    mat_angle_product(mat_angle_product<0) = NaN;   % 90°
    mat_angle_product(mat_angle_product>=0)=1;      % 180°
    
    % combine angle_sum and angle_product
    mat_angle = mat_angle_sum.*mat_angle_product;
    
    % delete data, if angle is too large or to small
    mat_angle_select = mat_angle;
    mat_angle_select(abs(mat_angle)<winkel_lower_bound)=NaN;
    mat_angle_select(abs(mat_angle)>winkel_upper_bound)=NaN;
    
    % define possible angles to "1"
    mat_poss_dam = mat_angle_select;
    mat_poss_dam(isnan(mat_poss_dam)==0)=1;
    
    
    %%% combine selected angles with selected distances
    
    mat_dist_dam = mat_dist_dam_all.*mat_poss_dam;
    
    
    if debug_on == 1
        figure('color',[1,1,1],'position',[1,1,2000,1000])
        subplot(2,3,1)
        img = imagesc(mat_angle_sum);
        set(img,'AlphaData',~isnan(mat_angle_sum))
        colorbar
        axis equal tight
        title('mat_angle_sum','interpreter','none')
        subplot(2,3,2)
        img = imagesc(mat_angle_product);
        set(img,'AlphaData',~isnan(mat_angle_product))
        colorbar
        axis equal tight
        title('mat_angle_product','interpreter','none')
        subplot(2,3,3)
        img = imagesc(mat_angle);
        set(img,'AlphaData',~isnan(mat_angle))
        colorbar
        axis equal tight
        title('mat_angle','interpreter','none')
        subplot(2,3,4)
        img = imagesc(mat_angle_select);
        set(img,'AlphaData',~isnan(mat_angle_select))
        colorbar
        axis equal tight
        title('mat_angle_select','interpreter','none')
        subplot(2,3,5)
        img = imagesc(mat_poss_dam);
        set(img,'AlphaData',~isnan(mat_poss_dam))
        colorbar
        axis equal tight
        title('mat_poss_dam','interpreter','none')
        subplot(2,3,6)
        img = imagesc(mat_dist_dam);
        set(img,'AlphaData',~isnan(mat_dist_dam))
        colorbar
        axis equal tight
        title('mat_dist_dam','interpreter','none')
    end
    
    %% save selected dam endpoints in matrix (refDAM)
    
    % determine location of shortest dam
    [row,col] = find(mat_dist_dam==min(min(mat_dist_dam)));
    
    % if there is a dam location
    if isempty(row)==0
        % index of the coordinates in the matrix
        idx = sub2ind(size(mat_dist_dam),row(1),col(1));
        
        % point location in the array
        point1 = array_raster_data_loc2(row(1),1);
        point2 = array_raster_data_loc1(col(1),1);

        
        %%% results
        
        % ids of the dam endpoints
        points = [point1,point2];
        
        % coordinates of the dam endpoints
        [points_coords(2,:),points_coords(1,:)] = ind2sub(size(refDAM_dem),[points(1);points(2)]);
        
        % points at the boundary of the raster behind point1 and point2
        points_boundary_coords = [matX2(point1),matX2(point2);...
            matY2(point1),matY2(point2)];
        
        % length of the minimum dam distance
        min_dist_dam = mat_dist_dam(idx(1,1));
        
        if debug_on == 1
            figure('color',[1,1,1],'position',[1,1,1000,1000])
            hold on
            img = imagesc(refDAM_waterline);
            set(img,'AlphaData',~isnan(refDAM_waterline))
            p1 = plot(points_coords(1,:),points_coords(2,:),'or');
            p2 = plot(points_boundary_coords(1,:),points_boundary_coords(2,:),'--ob');
            p3 = plot(point_xy(1),point_xy(2),'*r');
            legend([p3,p1,p2],'river point','dam end-points','points at the boundary behind the dam end-points','location','best')
            colorbar
            axis equal tight
            title('dam point positions on the shore line','interpreter','none')
            grid on
            box on
            keyboard
        end
        
        % stop loop (and exit function), if a dam was found (if not, a dam with lower
        % height is tried)
        if river_points(i_point).exit_code ~= 5
            if dam_height == dam_height_max
                % dam possible with requested height
                river_points(i_point).exit_code = 0;
                river_points_all(i_point_all).exit_code = 0;
            else
                % dam possible, but height is lower than the requested one
                river_points(i_point).exit_code = -1;
                river_points_all(i_point_all).exit_code = -1;
            end
        end
        river_points(i_point).dam_height = dam_height;
        river_points(i_point).dam_top_elev = river_points(i_point).dem + dam_height;
        river_points(i_point).dam_length = min_dist_dam;
        river_points(i_point).refDAM_LLcorner = [Xmin;Ymin];    % [X ; Y]
        river_points(i_point).refDAM_size = refDAM_size;
        
        river_points(i_point).refDAM_damEndpoints_coords = points_coords;   % [X ; Y]
        river_points(i_point).refDAM_damBoundarypoints_coords = points_boundary_coords;     % [X ; Y]
        river_points(i_point).refGLOBAL_damEndpoints_coords = points_coords + ...
            river_points(i_point).refDAM_LLcorner;    % [X ; Y]
        river_points(i_point).refGLOBAL_damBoundarypoints_coords = points_boundary_coords + ...
            river_points(i_point).refDAM_LLcorner;    % [X ; Y]
        
        river_points_all(i_point_all).refGLOBAL_damEndpoints_coords = river_points(i_point).refGLOBAL_damEndpoints_coords;
        
        
        % check if dam meets the criteria
        dam_length = sqrt((river_points(i_point).refDAM_damEndpoints_coords(1,2)-river_points(i_point).refDAM_damEndpoints_coords(1,1))^2 +...
            (river_points(i_point).refDAM_damEndpoints_coords(2,2)-river_points(i_point).refDAM_damEndpoints_coords(2,1))^2)*grids.cellsize;    % in m
                
        if river_points(i_point).exit_code<=0
            % exit function if valid dam was found
            break
        else
            if dam_height == dam_height_min
                river_points(i_point).exit_code = 6;  % no dam of this height was found for the point
                river_points_all(i_point_all).exit_code = 6;  % no dam of this height was found for the point
                return
            end
        end
    else
        if abs(dam_height-dam_height_min) <= dam_dist_eval
            river_points(i_point).exit_code = 6;  % no dam of this height was found for the point
            river_points_all(i_point_all).exit_code = 6;  % no dam of this height was found for the point
            return
        end
    end
    
end


%% exclude neighbors (if defined by user)
if river_points(i_point).exit_code<=0
river_points_all_id = [river_points_all.id];
river_points_id = [river_points.id];
neighbor_downstream = river_points_all(i_point_all).nachfolger;
    neighbor_upstream = river_points_all(i_point_all).vorgaenger;
    for n = 1:floor(neighbors_exclude_distance/grids.cellsize)
        if isempty(neighbor_downstream)==0 && isnan(neighbor_downstream)==0
            if river_points_all(river_points_all_id==neighbor_downstream).exit_code <= 0
                river_points_all(river_points_all_id==neighbor_downstream).exit_code = 5;
                i_river_point = find(river_points_id==neighbor_downstream);
                if isempty(i_river_point)==0
                    river_points(i_river_point).exit_code = 5;
                end
            end
            neighbor_downstream = river_points_all(river_points_all_id==neighbor_downstream).nachfolger;
        end
        if isempty(neighbor_upstream)==0
            if river_points_all(river_points_all_id==neighbor_upstream(1)).exit_code <= 0
                river_points_all(river_points_all_id==neighbor_upstream(1)).exit_code = 5;
                i_river_point = find(river_points_id==neighbor_upstream(1));
                if isempty(i_river_point)==0
                    river_points(i_river_point).exit_code = 5;
                end
            end
            neighbor_upstream = river_points_all(river_points_all_id==neighbor_upstream(1)).vorgaenger;
        end   
    end
end
