function fun_9n2_save_basins_characteristics_as_xlsx
% ######  save basin-curves  ######
% function to save curves (depth-storage-area) of the basins in a
% Excel-file (can be opened without MATLAB)
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% August 2019; Last revision: 22-Dez-2019

global grids basins_selected
warning off

% define filename
filename = 'basin_characteristics.xlsx';
delete(filename)

% define parameters
row_names = {'catchment area','dam height','dam axis length',...
    'dam volume','basin volume','specific volume',...
    'basin area','non-suited area'};

% combine data in matrix
for i_dam = 1:length(basins_selected)
i_selected = find(basins_selected(i_dam).curve_potential_dam_basin(:,1)==basins_selected(i_dam).dam_height);
col_names{i_dam} = sprintf('basin %d',i_dam);
data_matrix(:,i_dam) = [basins_selected(i_dam).curve_potential_dam_basin(i_selected,15),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,1),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,3),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,2),...    
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,4),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,6),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,5),...
    basins_selected(i_dam).curve_potential_dam_basin(i_selected,14)*(grids.cellsize)^2]';

end

% write general data
xlswrite(filename,row_names','general','A2');
xlswrite(filename,col_names,'general','B1');
xlswrite(filename,data_matrix,'general','B2');

% define names and data for curves (one sheet per basin)
col_names_curves = {'basin volume [m³]','water depth [m]','basin area [m²]'};
for i_dam = 1:length(basins_selected)
    xlswrite(filename,col_names_curves,sprintf('basin%d',i_dam),'A1');
    xlswrite(filename,basins_selected(i_dam).ShA_dam,sprintf('basin%d',i_dam),'A2');
end