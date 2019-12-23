function fun_4_determine_dam_characteristics(i_point)
% ######  determine dam characteristics  ######
% function to determine the geometry of the dam (volume, height, ...)
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global river_points debug_on grids dam_crest_width dam_slope_m
global dam_dist_eval dam_height_min dam_heights_eval
global point_boundaries

%% define and select input data

refDAM_dem = grids.dem(river_points(i_point).refDAM_LLcorner(2):...
    river_points(i_point).refDAM_LLcorner(2)+river_points(i_point).refDAM_size(2)-1,...
    river_points(i_point).refDAM_LLcorner(1):...
    river_points(i_point).refDAM_LLcorner(1)+river_points(i_point).refDAM_size(1)-1);

pkt_damcoord = river_points(i_point).refDAM_damBoundarypoints_coords;

% water depth for a dam with maximum height
refDAM_waterdepth_max = refDAM_dem - (river_points(i_point).dem + river_points(i_point).dam_height);


%% determine X, Y and elevation information for the dam
%
% --> interpolate points on the dam (dam endpoints and river point)
%

% separate analysis of the dam segments
% (left endpoint to river point & right endpoint to river point)
for dam_segment = 1:2
    
    % length of the dam (in x or y direction, no diagonal)
    dam_length{dam_segment} = max(abs(pkt_damcoord(1,1)-pkt_damcoord(1,2)),abs(pkt_damcoord(2,1)-pkt_damcoord(2,2)));
    
    % interpolate x-coordinates (column numbers)
    if pkt_damcoord(1,1)==pkt_damcoord(1,2) % if the dam is absolutely vertical
        dampoints_segment_X{dam_segment} = repmat(pkt_damcoord(1,1),1,dam_length{dam_segment}+1);
    else
        dampoints_segment_X{dam_segment} = round(pkt_damcoord(1,1):(pkt_damcoord(1,2)-pkt_damcoord(1,1))/dam_length{dam_segment}:pkt_damcoord(1,2));
    end
    
    % interpolate y-coordinates (row numbers)
    if pkt_damcoord(2,1)==pkt_damcoord(2,2) % if the dam is absolutely horizontal
        dampoints_segment_Y{dam_segment} = repmat(pkt_damcoord(2,1),1,dam_length{dam_segment}+1);
    else
        dampoints_segment_Y{dam_segment} = round(pkt_damcoord(2,1):(pkt_damcoord(2,2)-pkt_damcoord(2,1))/dam_length{dam_segment}:pkt_damcoord(2,2));
    end
    
end

% combine both dam segments
dampoints_X = [dampoints_segment_X{1},dampoints_segment_X{2}];
dampoints_Y = [dampoints_segment_Y{1},dampoints_segment_Y{2}];

dam_axis = unique([dampoints_X;dampoints_Y]','rows','stable');

% sort dam axis coordinates
if (dam_axis(1,1) > dam_axis(end,1) && dam_axis(1,2) > dam_axis(end,2)) ||...
        (dam_axis(1,1) < dam_axis(end,1) && dam_axis(1,2) < dam_axis(end,2))
    % positive axis slope of the dam
    dam_axis = sortrows(dam_axis,2);
    dam_axis = sortrows(dam_axis,1);
else 
    % negative axis slope of the dam
    dam_axis = sortrows(dam_axis,-2);
    dam_axis = sortrows(dam_axis,1);
end
dam_axis = dam_axis';

% index of the dam coordinates in the GRID
index_dampoints = sub2ind(size(refDAM_waterdepth_max),dam_axis(2,:),dam_axis(1,:));

% get relative height of the dam points (difference to the top elevation of the dam)
dampoints_height = refDAM_waterdepth_max(index_dampoints);

% combine X, Y and elevation in one variable
dam_axis(3,:) = dampoints_height;



%% determine cross-section of the dam

dam_heights_eval = [river_points(i_point).dam_height:-dam_dist_eval:dam_height_min];
dam_elevation_eval = dam_heights_eval+river_points(i_point).dem;

i_top = 0;
i_seite = 0;

% angle of the dam with respect to the horizontal line
beta = abs(atan((river_points(i_point).refDAM_damEndpoints_coords(2,1)-river_points(i_point).refDAM_damEndpoints_coords(2,2))/...
            (river_points(i_point).refDAM_damEndpoints_coords(1,1)-river_points(i_point).refDAM_damEndpoints_coords(1,2))));
        
%%% widening of the dam (crest and slope) according to user information
% definition of the widening range (total width of the dam in flow direction)
dam_crest_multiplication = round(dam_crest_width/grids.cellsize) + round(river_points(i_point).dam_height/grids.cellsize)*dam_slope_m*4;
dam_crest_Y_temp2 = [];
dam_crest_X_temp2 = [];
dam_crest_Z_temp2 = [];
for steps = (-1)*floor(dam_crest_multiplication/2):floor(dam_crest_multiplication/2)

    dam_crest_X_temp2 = [dam_crest_X_temp2, dam_axis(1,:)+steps, dam_axis(1,:)];
    dam_crest_Y_temp2 = [dam_crest_Y_temp2, dam_axis(2,:), dam_axis(2,:)+steps];
    
    if abs(steps) <= floor(dam_crest_width/grids.cellsize/2)
        elev_dam_x = dam_elevation_eval';
        elev_dam_y = dam_elevation_eval';
        i_top = i_top+1;
    else
        Lx = sin(beta)*(abs(steps*grids.cellsize)-floor(dam_crest_width/grids.cellsize/2))/dam_slope_m;
        Ly = cos(beta)*(abs(steps*grids.cellsize)-floor(dam_crest_width/grids.cellsize/2))/dam_slope_m;
        elev_dam_x = dam_elevation_eval' - Lx;
        elev_dam_y = dam_elevation_eval' - Ly;
        i_seite = i_seite+1;
    end
    dam_crest_Z_temp2 = [dam_crest_Z_temp2, repmat(elev_dam_x,1,length(dam_axis(1,:))),repmat(elev_dam_y,1,length(dam_axis(2,:)))];
    if steps == 0
        dam_crest_axis(1,:) = [dam_axis(1,:)];
        dam_crest_axis(2,:) = [dam_axis(2,:)];
        dam_crest_axis(3:2+length(dam_elevation_eval),:) = ...
            repmat(dam_elevation_eval',1,length(dam_axis(1,:)));
    end
end
dam_crest_temp = [dam_crest_X_temp2;dam_crest_Y_temp2;dam_crest_Z_temp2]';

% exclude double coordinates in dam crest
dam_crest = unique(dam_crest_temp(:,1:2),'rows');
dam_crest_temp = sortrows(dam_crest_temp,-3);
[Lia,Locb] = ismember(dam_crest(:,1:2),dam_crest_temp(:,1:2),'rows');
dam_crest(Lia,3:size(dam_crest_temp,2)) = dam_crest_temp(Locb,3:end);


%%% more precise way for the description (but larger computational effort)
% [dam_crest,~,k] = unique(dam_crest_temp(:,1:2),'rows');
% for nu = 1:numel(dam_crest)
%     indexToThisUniqueValue = (nu==k)';
%     dam_crest(nu,3:size(dam_crest_temp,2)) = mean(dam_crest_temp(indexToThisUniqueValue,3:end),1);
% end

dam_crest = dam_crest';

% plausibility check (dam within raster)
dam_crest = dam_crest(:,dam_crest(2,:) > max(1,min(river_points(i_point).refDAM_damEndpoints_coords(2,:))-point_boundaries/10) ...
    & dam_crest(1,:) > max(1,min(river_points(i_point).refDAM_damEndpoints_coords(1,:))-point_boundaries/10) ...
    & dam_crest(2,:) < min(size(refDAM_waterdepth_max,1),max(river_points(i_point).refDAM_damEndpoints_coords(2,:))+point_boundaries/10) ...
    & dam_crest(1,:) < min(size(refDAM_waterdepth_max,2),max(river_points(i_point).refDAM_damEndpoints_coords(1,:))+point_boundaries/10));


%% determine dam heights, elevations and lengths

% dam crest
index_dam_surface = sub2ind(size(refDAM_dem),dam_crest(2,:),dam_crest(1,:));

dam_heights = [dam_crest(1,:);dam_crest(2,:);
    dam_crest(3:end,:)-refDAM_dem(index_dam_surface)];
dam_heights = dam_heights(:,dam_heights(3,:)>=0);
refDAM_dam_coords = dam_heights(1:2,:);
dam_heights = dam_heights(3:end,:);

index_dam_surface = sub2ind(size(refDAM_dem),refDAM_dam_coords(2,:),refDAM_dam_coords(1,:));
dam_heights(dam_heights<0) = NaN;
dam_elevations = dam_heights+refDAM_dem(index_dam_surface);

% dam axis
index_dam_crest_axis = sub2ind(size(refDAM_dem),dam_crest_axis(2,:),dam_crest_axis(1,:));
dam_axis_heights = [nan(1,length(index_dam_crest_axis));index_dam_crest_axis;dam_crest_axis(1,:);dam_crest_axis(2,:);
    dam_crest_axis(3:end,:)-refDAM_dem(index_dam_crest_axis)]';
% dam_axis_heights = sortrows(dam_axis_heights,2);
dam_axis_heights = dam_axis_heights';
dam_axis_heights = dam_axis_heights(:,dam_axis_heights(5,:)>=0);

% determine the lengths of the dam axis segments (to plot dam crosssection)
alpha_dam = atan([(dam_axis_heights(4,end)-dam_axis_heights(4,1))./(dam_axis_heights(3,end)-dam_axis_heights(3,1))]);
alpha_segment = [alpha_dam,atan([(dam_axis_heights(4,2:end)-dam_axis_heights(4,1:end-1))...
    ./(dam_axis_heights(3,2:end)-dam_axis_heights(3,1:end-1))])];
alpha = (alpha_dam-alpha_segment);
distances_segments = [0,((dam_axis_heights(3,2:end)-dam_axis_heights(3,1:end-1)).^2 + ...
    (dam_axis_heights(4,2:end)-dam_axis_heights(4,1:end-1)).^2).^(1/2)]*grids.cellsize;

dam_axis_heights(1,:) = cos(alpha).*distances_segments;

refDAM_dam_axis_coords = dam_axis_heights(3:4,:);
dam_segment_length = repmat((dam_axis_heights(1,:)),size(dam_axis_heights,1)-4,1);
dam_segment_length(dam_axis_heights(5:end,:)<=0) = 0;
dam_axis_lengths = sum(dam_segment_length,2);
dam_axis_heights = dam_axis_heights(5:end,:);


%% determine dam volumes

dam_volumes(:,1) = sum(dam_heights.*(grids.cellsize)^2,2,'omitnan');


%% save variables

river_points(i_point).refDAM_dam_coords = refDAM_dam_coords;
river_points(i_point).dam_heights_eval = dam_heights_eval;%(lines);
river_points(i_point).dam_heights = dam_heights;
river_points(i_point).dam_elevations = dam_elevations;
river_points(i_point).dam_volumes = dam_volumes;

river_points(i_point).refDAM_dam_axis_coords = refDAM_dam_axis_coords;
river_points(i_point).dam_axis_lengths = dam_axis_lengths;
river_points(i_point).dam_axis_lengths_segments = dam_segment_length(1,:);
river_points(i_point).dam_axis_heights = dam_axis_heights;


if debug_on == 1
    
    %% top view of dam
    
    i_dam = 1;    
    
    % definition of colormap
    anzahl = 10;
    orange = [177, 108, 37]/255;
    rot = [141, 36, 37]/255;
    blau1 = [105, 151, 201]/255;
    blau4 = [36, 68, 101]/255;
    orange_rot = [linspace(orange(1),rot(1),anzahl)',linspace(orange(2),rot(2),anzahl)',linspace(orange(3),rot(3),anzahl)'];
    blau14 = [linspace(blau1(1),blau4(1),anzahl)',linspace(blau1(2),blau4(2),anzahl)',linspace(blau1(3),blau4(3),anzahl)'];
    cmap = [flipud(blau14);orange_rot];
    % definition of the range of the colormap (based on local water depth grid)
    max_abs_val = max(max(abs(refDAM_waterdepth_max)));
    
    % FIGURE
    figure('color',[1,1,1],'position',[1,1,1000,1000])

    % plot grid with colorbar
    ax(1) = axes;
    hold on
    imagesc(refDAM_waterdepth_max)
    caxis([-max_abs_val,max_abs_val])
    
    % plot points of dam and river points in the basin
    ax(2) = axes;
    hold on
    p2 = plot(river_points(i_point).refDAM_damEndpoints_coords(1,:),river_points(i_point).refDAM_damEndpoints_coords(2,:),'or');
    p3 = plot(river_points(i_point).refDAM_damBoundarypoints_coords(1,:),river_points(i_point).refDAM_damBoundarypoints_coords(2,:),'og');
    p4 = scatter(river_points(i_point).refDAM_dam_coords(1,river_points(i_point).dam_heights(i_dam,:)>0),...
        river_points(i_point).refDAM_dam_coords(2,river_points(i_point).dam_heights(i_dam,:)>0),15,...
        river_points(i_point).dam_heights(i_dam,river_points(i_point).dam_heights(i_dam,:)>0),'filled');
    p5 = scatter(river_points(i_point).refDAM_dam_axis_coords(1,:),river_points(i_point).refDAM_dam_axis_coords(2,:),'.k');
    p1 = plot(river_points(i_point).col-river_points(i_point).refDAM_LLcorner(1)+1,...
        river_points(i_point).row-river_points(i_point).refDAM_LLcorner(2)+1,'*r');
    set(ax(2),'color','none')
    % legend
    legend([p1,p2,p3,p4,p5],'river point of the dam','dam endpoints','points at the refDAM-boundaries','dam points','dam axis','fontsize',11,'location','best')
    
    % Give each one its own colormap
    colormap(ax(1),cmap)
    colormap(ax(2),'parula')
    % add colorbars
    cb1 = colorbar(ax(1),'location','eastoutside');
    ylabel(cb1,'elevation relative to the elevation of the dam top (= water depth) (m)','fontsize',11)
    cb2 = colorbar(ax(2),'location','northoutside');
    ylabel(cb2,'dam height (m)','fontsize',11)
    
    for i = 1:2
        if i>1
            ax(i).Visible = 'off';
            ax(i).XTick = [];
            ax(i).YTick = [];
        end
        set(ax(i),'position',[.1,.1,.8,.8],'xlim',[0,size(refDAM_waterdepth_max,2)],'ylim',[0,size(refDAM_waterdepth_max,1)],'fontsize',11)
    end
    linkaxes(ax(1),ax(2),'xy')
    
    
    %% dam cross section
    i_dam = 1;
    
    figure('color',[1,1,1])
    plot([0,cumsum(dam_segment_length(i_dam,dam_segment_length(i_dam,:)>0)),max(cumsum(dam_segment_length(i_dam,dam_segment_length(i_dam,:)>0)))],...
        [0,(-1)*dam_axis_heights(i_dam,dam_segment_length(i_dam,:)>0),0])
    xlabel('dam length (m)')
    ylabel('dam height (m)')
    legend('dam cross section','location','best')
    set(gca,'xlim',[0,max(cumsum(dam_segment_length(i_dam,dam_segment_length(i_dam,:)>0)))])
    grid on
    grid minor
    
    %% dam volume
    
    [~,col] = find(dam_axis_heights==max(max(dam_axis_heights)));
    heights = dam_axis_heights(:,col);
        
    figure('color',[1,1,1])
    plot(heights,dam_volumes)
    xlabel('maximum dam height (m)')
    ylabel('dam volume (m³)')
    grid on
    grid minor
    
    keyboard
end