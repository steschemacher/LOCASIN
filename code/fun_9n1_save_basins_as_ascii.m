function fun_9n1_save_basins_as_ascii
% ######  save spatial distribution of basins  ######
% function to save the basin characteristics as ascii-file which can be
% imported to GIS-systems. The dam heights are defined positive, while
% water depths are defined negative
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global basins_selected


for i_dam = 1:length(basins_selected)
    
    BASIN_grid = basins_selected(i_dam).refBASIN_depth_dam_crosssection * (-1);
        
    % convert dam coordinates to raster
    dam_x = basins_selected(i_dam).refGK_BASIN_x;
    dam_y = basins_selected(i_dam).refGK_BASIN_y;
    grid.cellsize = dam_x(2)-dam_x(1);
    DAM_grid = nan(length(dam_y),length(dam_x));
    index_dam = sub2ind(size(DAM_grid),(max(dam_y)-min(dam_y))/grid.cellsize - (basins_selected(i_dam).refGK_dam_coords(2,:)-min(dam_y))/grid.cellsize+1,...
        (basins_selected(i_dam).refGK_dam_coords(1,:)-min(dam_x))/grid.cellsize+1);
    DAM_grid(index_dam) = basins_selected(i_dam).dam_heights(1,:);
    
    grid.data = DAM_grid;
    grid.data(isnan(BASIN_grid)==0) = BASIN_grid(isnan(BASIN_grid)==0);
    grid.data = grid.data;
    grid.ncols = length(dam_x);
    grid.nrows = length(dam_y);
    grid.xll = min(dam_x);
    grid.yll = min(dam_y);
    
    grid.data(isnan(grid.data)==1) = -9999;
    
    % write the data in a txt-file
    file_txt = sprintf('basin_%03d.txt',i_dam);

    txt_file = fopen(file_txt,'w');
    fprintf(txt_file,sprintf('ncols         %d \n',grid.ncols));
    fprintf(txt_file,sprintf('nrows         %d \n',grid.nrows));
    fprintf(txt_file,sprintf('xllcorner     %d \n',grid.xll));
    fprintf(txt_file,sprintf('yllcorner     %d \n',grid.yll));
    fprintf(txt_file,sprintf('cellsize      %d \n',grid.cellsize));
    fprintf(txt_file,sprintf('NODATA_value  -9999 \n'));
    fclose(txt_file);
    dlmwrite(file_txt,grid.data,'delimiter','	','precision',7,'-append');
end
