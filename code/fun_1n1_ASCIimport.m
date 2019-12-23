function [grid] = fun_1n1_ASCIimport(pfad,name,x_range,y_range)
% ######  import ASCII data  ######
% import ASCII data into matlab including the coordinates of the reference
% system and the size information
%
% Inputs:
%    pfad: directory at which the data is located
%    name: name of data which needs to be loaded
%
% Outputs: grid struct
%    ncols: number of columns of the grid
%    nrows: number of rows of the grid
%    xll:   x-coordinate of the lower left corner (in meters)
%    yll:   y-coordinate of the lower left corner (in meters)
%    cellsize: cellsize in meters
%    x:     x-coordinates of every column
%    y:     y-coordinates of every row
%    X:     x-coordinates of the whole matrix (size: ncols*nrows)
%    Y:     y-coordinates of the whole matrix (size: ncols*nrows)
%
% Other m-files required: -
% functions: -
% MAT-files required: -
%
% Author: Sonja Tesschemacher
% email: sonja.teschemacher@tum.de
% August 2016; Last revision: 26-Oct-2019


%% read data
cd(pfad)

if isempty(x_range)==0
    %%% read only part of input data
    content = fileread( name ) ;
    headers = textscan(content,'%s',12);
    % read header data
    ncols = str2double(headers{1}(2));
    nrows = str2double(headers{1}(4));
    xll = str2double(headers{1}(6));
    yll = str2double(headers{1}(8));
    cellsize = str2double(headers{1}(10));
    nodata = str2double(headers{1}(12));
    
    % define rows and columns
    x_min = x_range(1);
    x_max = x_range(2);
    y_min = y_range(1);
    y_max = y_range(2);
    
    headerlines = 6;
    rows = [headerlines + nrows - round((y_min-yll)/cellsize) , headerlines + nrows - round((y_max-yll)/cellsize)];
    columns = [round((x_min-xll)/cellsize) , round((x_max-xll)/cellsize)];
    
    % select section
    data = textscan( content,[repmat('%*f', 1, columns(1)),repmat('%f',1,columns(2)-columns(1)+1),'%*[^\n]'],...
        rows(1)-rows(2)+1,'HeaderLines', rows(2)) ;
    data = cat(2,data{:});
    
    
    %%% reorder data of the header
    
    % size of data
    grid.ncols = size(data,2);
    grid.nrows = size(data,1);
    
    % coordinate reference
    grid.xll = x_min;
    grid.yll = y_min;
    
    % cellsize
    grid.cellsize = cellsize;
    
    % arrays and matrices to plot with reference system
    grid.x = grid.xll:grid.cellsize:grid.xll+grid.cellsize*(grid.ncols-1);
    grid.y = grid.yll+grid.cellsize*(grid.nrows-1):-grid.cellsize:grid.yll;
    [grid.X,grid.Y] = meshgrid(grid.x,grid.y);
        
    %%% include data
    grid.data = data;
    grid.data(grid.data==-9999) = NaN;
    
else
    %%% read complete files
    data_raw = importdata(name,' ',6);  % data separated by blanks
    
    if iscell(data_raw)==1  % data separated by tabs
        data_raw = importdata(name,'	',6);
    end
    
    
    %%% reorder data of the header
    
    % size of data
    str_temp = char(data_raw.textdata(1,1));
    grid.ncols = str2num(str_temp(6:length(str_temp)));
    clear str_temp
    str_temp = char(data_raw.textdata(2,1));
    grid.nrows = str2num(str_temp(6:length(str_temp)));
    clear str_temp
    
    % coordinate reference
    str_temp = char(data_raw.textdata(3,1));
    grid.xll = str2num(str_temp(10:length(str_temp)));
    clear str_temp
    str_temp = char(data_raw.textdata(4,1));
    grid.yll = str2num(str_temp(10:length(str_temp)));
    clear str_temp
    
    % cellsize
    str_temp = char(data_raw.textdata(5,1));
    grid.cellsize = str2num(str_temp(9:length(str_temp)));
    clear str_temp
    
    % arrays and matrices to plot with reference system
    grid.x = grid.xll:grid.cellsize:grid.xll+grid.cellsize*(grid.ncols-1);
    grid.y = grid.yll+grid.cellsize*(grid.nrows-1):-grid.cellsize:grid.yll;
    [grid.X,grid.Y] = meshgrid(grid.x,grid.y);
        
    %%% include data
    grid.data = data_raw.data;
    grid.data(grid.data==-9999) = NaN;
    
end


