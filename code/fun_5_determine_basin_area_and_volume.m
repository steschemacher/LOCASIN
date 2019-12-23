function fun_5_determine_basin_area_and_volume(i_point,i_point_all)
% ######  determine basin characteristics  ######
% function to determine the inundation area of the basin and the subsequent
% depth distribution and volume
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 16-Dez-2019

global river_points grids dam_crest_width debug_on
global river_points_all dam_length_max dam_height_buffer


% define boundary coordinates
Xmin = river_points(i_point).refGLOBAL_basin_boundaries(1,1);
Xmax = river_points(i_point).refGLOBAL_basin_boundaries(1,2);
Ymin = river_points(i_point).refGLOBAL_basin_boundaries(2,1);
Ymax = river_points(i_point).refGLOBAL_basin_boundaries(2,2);

n_tries_max = 4;
n_tries = 0;
shift_border = round(dam_length_max/grids.cellsize/n_tries_max);
dam_heights_eval = river_points(i_point).dam_heights_eval;
dam_elevations_eval = dam_heights_eval + river_points(i_point).dem;
refGLOBAL_dam_coords = river_points(i_point).refDAM_dam_coords + (river_points(i_point).refDAM_LLcorner - [1;1]);
dam_heights = river_points(i_point).dam_heights;


% determine point in the upstream area of the dam
% (to calculate the flooding area)
n_points_upstream = min(ceil(dam_crest_width/grids.cellsize)*2,length(river_points(i_point).basin_points));


if n_points_upstream~=0
    point_upstream = river_points(i_point).basin_points(n_points_upstream);    % point id

    while 1
        
        % select DEM        
        refBASIN_DEM = grids.dem(Ymin:Ymax,Xmin:Xmax);
        refBASIN_LLcorner = [Xmin;Ymin];
        [size_row,size_col] = size(refBASIN_DEM);
        numel_raster = numel(refBASIN_DEM);
        
        
        %%% generate 3D-Raster with DAM surface and DEM
        % dam points in refBASIN konvertieren
        refBASIN_dam_coords = refGLOBAL_dam_coords - (refBASIN_LLcorner -[1;1]);
        number_heights = length(dam_heights_eval);
        refBASIN_DAM = zeros(size_row,size_col,number_heights);
        index_dam = sub2ind([size_row,size_col],refBASIN_dam_coords(2,:),refBASIN_dam_coords(1,:));
        index_dam_mat = reshape((index_dam + [([1:number_heights]-1)*(numel_raster)]')',1,[]);
        refBASIN_DAM(index_dam_mat) = reshape(dam_heights',1,[]);
        refBASIN_DEM = repmat(refBASIN_DEM,1,1,number_heights);
        refBASIN_DAM(isnan(refBASIN_DAM)==1) = 0;
        refBASIN_DEM_DAM = refBASIN_DAM + refBASIN_DEM;
        
        refBASIN_waterlevel_elevation = reshape(repmat(dam_elevations_eval',1,numel_raster)',size_row,size_col,[]);
        
        refBASIN_water_depths = refBASIN_DEM_DAM-refBASIN_waterlevel_elevation;
        LAND = refBASIN_water_depths(:,:,1)+dam_height_buffer;
        LAND(LAND<=0) = 0;
        LAND(LAND>0) = 1;
        WASSER = refBASIN_water_depths+dam_height_buffer;
        WASSER(WASSER>=0) = NaN;
        WASSER(WASSER<0) = 1;
        
        %%% determine lake area by extending around upstream point
        % determine the coordinates of the point (in the refBASIN reference
        % system)
        punkt = [river_points_all([river_points_all.id]==point_upstream).col;...
            river_points_all([river_points_all.id]==point_upstream).row]-(refBASIN_LLcorner-[1;1]);
       
        while LAND(punkt(2),punkt(1))==1
            n_points_upstream = n_points_upstream + 2;
            if length(river_points(i_point).basin_points)>=n_points_upstream
                point_upstream = river_points(i_point).basin_points(n_points_upstream);
                punkt = [river_points_all([river_points_all.id]==point_upstream).col;...
                    river_points_all([river_points_all.id]==point_upstream).row]-(refBASIN_LLcorner-[1;1]);
            else
                % only river points on the shoreline or dam could be found
                river_points(i_point).exit_code = 9;  
                river_points_all(i_point_all).exit_code = 9;
                return
            end
        end
        BASIN_max = encodem(LAND,[punkt(2),punkt(1),3]);
        BASIN_max(BASIN_max~=3) = NaN;
        BASIN_max(BASIN_max==3) = 1;
        
        BASIN = repmat(BASIN_max,1,1,number_heights);
        BASIN = BASIN.*WASSER;
        
        clear LAND WASSER refBASIN_water_depths
        
        % test if basin reaches the borders of the grid
        BORDER = zeros(size(BASIN_max));
        BORDER(:,1) = 1;    % left border
        BORDER(:,end) = 2;  % right border
        BORDER(1,:) = 3;    % top border (xmin)
        BORDER(end,:) = 4;  % bottom border (xmax)
        n_tries = n_tries+1;
        
        if sum(BASIN_max(BORDER>0),'omitnan')==0
            break
        else
            if n_tries == n_tries_max
                if sum(BASIN_max(BORDER>0),'omitnan')~=0
                    for i_test_height = 1:length(dam_heights_eval)
                        BASIN_max = BASIN(:,:,i_test_height);
                        if sum(BASIN_max(BORDER>0),'omitnan')==0
                            dam_heights_eval = dam_heights_eval(i_test_height:end);
                            refBASIN_DEM_DAM = refBASIN_DEM_DAM(:,:,i_test_height:end);
                            refBASIN_DEM = refBASIN_DEM(:,:,i_test_height:end);
                            BASIN = BASIN(:,:,i_test_height:end);
                            refBASIN_DAM = refBASIN_DAM(:,:,i_test_height:end);
                            refBASIN_waterlevel_elevation = refBASIN_waterlevel_elevation(:,:,i_test_height:end); 
                            river_points(i_point).dam_heights = dam_heights(i_test_height:end,:);
                            river_points(i_point).dam_elevations = river_points(i_point).dam_elevations(i_test_height:end,:);
                            river_points(i_point).dam_volumes = river_points(i_point).dam_volumes(i_test_height:end,1);
                            river_points(i_point).dam_axis_lengths = river_points(i_point).dam_axis_lengths(i_test_height:end,:);
                            river_points(i_point).dam_axis_heights = river_points(i_point).dam_axis_heights(i_test_height:end,:);
                            
                            % dam heights needed to be reduced, because the
                            % basin touched the borders
                            river_points(i_point).exit_code = -3;
                            river_points_all(i_point_all).exit_code = -3;
                            break
                        end
                        
                    end
                end
                if sum(BASIN_max(BORDER>0),'omitnan')~=0
                    river_points(i_point).exit_code = 10;  % basin touches the borders
                    river_points_all(i_point_all).exit_code = 10;  % basin touches the borders
                    return
                else
                    break
                end
            end
            
            if sum(BASIN_max(BORDER==1),'omitnan')>0
                Xmin = Xmin - shift_border;
                if Xmin < 1
                    Xmin = 1;
                    n_tries = n_tries_max-1; 
                end
            end
            if sum(BASIN_max(BORDER==2),'omitnan')>0
                Xmax = Xmax + shift_border;
                if Xmax > length(grids.x)
                    Xmax = length(grids.x);
                    n_tries = n_tries_max-1; 
                end
            end
            if sum(BASIN_max(BORDER==4),'omitnan')>0
                Ymax = Ymax + shift_border;
                if Ymax > length(grids.y)
                    Ymax = length(grids.y);
                    n_tries = n_tries_max-1; 
                end
            end
            if sum(BASIN_max(BORDER==3),'omitnan')>0
                Ymin = Ymin - shift_border;
                if Ymin < 1
                    Ymin = 1;
                    n_tries = n_tries_max-1; 
                end
            end
            
            
            
        end
    end
            
else
    
    % no upstream point available (required for calculation)
    river_points(i_point).exit_code = 8;
    river_points_all(i_point_all).exit_code = 8;
    return
end
    
        
    refBASIN_waterlevel_elevation(BASIN~=1) = NaN;
    
    refBASIN_depth_dam_crosssection = refBASIN_waterlevel_elevation - refBASIN_DEM_DAM;
    basin_volume(:,1) = squeeze(sum(sum(refBASIN_depth_dam_crosssection,1,'omitnan'),2,'omitnan'))*(grids.cellsize)^2;
    basin_area(:,1) = sum(sum(refBASIN_depth_dam_crosssection>=0,1,'omitnan'),2,'omitnan')*(grids.cellsize)^2;
    
    % basin depth and area information (dam is vertical line)
    refBASIN_depth_wall_crosssection = refBASIN_waterlevel_elevation - refBASIN_DEM;
    basin_volume(:,2) = squeeze(sum(sum(refBASIN_depth_wall_crosssection,1,'omitnan'),2,'omitnan'))*(grids.cellsize)^2;
    basin_area(:,2) = sum(sum(refBASIN_depth_wall_crosssection>=0,1,'omitnan'),2,'omitnan')*(grids.cellsize)^2;
    
    % combine final basin and dam for maximum water depth
    DAM_max = refBASIN_DAM(:,:,1);
    BASIN_max(isnan(BASIN_max)==1) = 0;
    BASIN_DAM_max = DAM_max+BASIN_max;
    BASIN_DAM_max(BASIN_DAM_max==0) = NaN;
    
    % clip raster files to the maximum lake extent
    sumX = mean(BASIN_DAM_max,1,'omitnan');
    sumY = mean(BASIN_DAM_max,2,'omitnan');
    Xmin = max(1,find(isnan(sumX)==0,1,'first')-1);
    Xmax = min(length(sumX),find(isnan(sumX)==0,1,'last')+1);
    Ymin = max(1,find(isnan(sumY)==0,1,'first')-1);
    Ymax = min(length(sumY),find(isnan(sumY)==0,1,'last')+1);
    
    refBASIN_depth_dam_crosssection = refBASIN_depth_dam_crosssection(Ymin:Ymax,Xmin:Xmax,:);
    refBASIN_depth_wall_crosssection = refBASIN_depth_wall_crosssection(Ymin:Ymax,Xmin:Xmax,:);
    refBASIN_LLcorner = refBASIN_LLcorner + [Xmin;Ymin] - [1;1];
    refGK_BASIN_x = grids.x(refBASIN_LLcorner(1):refBASIN_LLcorner(1)+size(refBASIN_depth_dam_crosssection,2)-1);
    refGK_BASIN_y = grids.y(refBASIN_LLcorner(2):refBASIN_LLcorner(2)+size(refBASIN_depth_dam_crosssection,1)-1);
    
    
    %% save results
    %     river_points(i_point).basin_points = [];
    river_points(i_point).dam_volumes(:,2) = basin_volume(:,2) - basin_volume(:,1);
    river_points(i_point).basin_volumes = basin_volume;
    river_points(i_point).basin_areas = basin_area;
    river_points(i_point).dam_heights_eval = dam_heights_eval;
    river_points(i_point).refBASIN_depth_dam_crosssection = refBASIN_depth_dam_crosssection;
    river_points(i_point).refBASIN_depth_wall_crosssection = refBASIN_depth_wall_crosssection;
    river_points(i_point).refGK_BASIN_x = refGK_BASIN_x;
    river_points(i_point).refGK_BASIN_y = refGK_BASIN_y;
    river_points(i_point).refBASIN_LLcorner = refBASIN_LLcorner;
    river_points(i_point).refGLOBAL_basin_boundaries = [];
   

    
    if debug_on == 1
        %% dam volumes, basin volumes and basin areas for different potential dam heights
        figure('color',[1,1,1])
        hold on
        plot(river_points(i_point).dam_volumes(:,1),river_points(i_point).dam_heights_eval,'.-')
        plot(river_points(i_point).basin_areas(:,1),river_points(i_point).dam_heights_eval,'.-')
        plot(river_points(i_point).basin_volumes(:,1),river_points(i_point).dam_heights_eval,'.-')
        ylabel('potential dam heights (m)')
        xlabel('volume (m³) or area (m²)')
        legend('dam volumes (m³)','basin areas (m³)','basin volumes (m³)','location','southeast')
        grid on
        grid minor
        
        %% digital elevation model with dam and lake
        anzahl = 10;
        anzahl2 = 10;
        orange = [177, 108, 37]/255;
        rot = [141, 36, 37]/255;
        blau1 = [105, 151, 201]/255;
        blau4 = [36, 68, 101]/255;
        orange_rot = [reshape(repmat(linspace(orange(1),rot(1),anzahl),anzahl2,1),[],1),...
            reshape(repmat(linspace(orange(2),rot(2),anzahl),anzahl2,1),[],1),...
            reshape(repmat(linspace(orange(3),rot(3),anzahl),anzahl2,1),[],1)];
        blau14 = [reshape(repmat(linspace(blau1(1),blau4(1),anzahl),anzahl2,1),[],1),...
            reshape(repmat(linspace(blau1(2),blau4(2),anzahl),anzahl2,1),[],1),...
            reshape(repmat(linspace(blau1(3),blau4(3),anzahl),anzahl2,1),[],1)];
        
        %%% DAM SELECTION %%%
        i_dam = 1;
        %%%%%%%%%%%%%%%%%%%%
        
        figure('color',[1,1,1],'position',[1,1,1000,1000])
        
        % plot DEM with dam
        ax(1) = axes;
        hold on
        img = imagesc(refBASIN_DEM_DAM(:,:,i_dam));
        caxis([mean(refBASIN_DEM_DAM(:,:,i_dam),'all')-8,mean(refBASIN_DEM_DAM(:,:,i_dam),'all')+2])
        axis equal
        
        % plot lake on top
        ax(2) = axes;
        hold on
        img = imagesc(Xmin:Xmax,Ymin:Ymax,refBASIN_depth_dam_crosssection(:,:,i_dam));
        set(img,'AlphaData',~isnan(refBASIN_depth_dam_crosssection(:,:,i_dam)))
        caxis([0,8])
        axis equal
        
        % change colormaps
        colormap(ax(1),orange_rot)
        colormap(ax(2),blau14)
        % add colorbars and labels
        cb1 = colorbar(ax(1),'location','eastoutside');
        ylabel(cb1,'ground elevation (m+NN)','fontsize',11)
        cb2 = colorbar(ax(2),'location','northoutside');
        ylabel(cb2,'water depth (m)','fontsize',11)
        
        % set axis limits
        for i = 1:2
            if i>1
                ax(i).Visible = 'off';
                ax(i).XTick = [];
                ax(i).YTick = [];
            end
            set(ax(i),'position',[.1,.1,.8,.8],'xlim',[0,size(refBASIN_DEM_DAM(:,:,i_dam),2)],'ylim',[0,size(refBASIN_DEM_DAM(:,:,i_dam),1)],'fontsize',11)
        end
        keyboard
    end
   