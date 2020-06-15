% ######  visual representation of the results  ######
% separate plotting of the figures from already calculated *.mat-results
%
% functions:    -
%
% Author: Sonja Teschemacher
% email: sonja.teschemacher@tum.de
% December 2019; Last revision: 15-Jun-2020

close all
global path_result
read_from_xlsx = 1;

if read_from_xlsx==1
    info_input = textread('define_input_directory_and_file_plotting.txt','%s');
    name_directory = info_input{1};
    try
        cd(name_directory)
    catch
        error(sprintf('The defined input directory is not existing:\n %s',name_directory))
    end
end
diary('logfile_plotting.txt')

try    
    %-- define status figure ---------------------------------------------%
    logo = imread('locasin_plotting.jpg');
    overview = figure('name','status','NumberTitle','off','visible','off');
    s = subplot(1,1,1);
    imshow(logo)
    size_position = get(overview,'Position');
    set(overview,'Position',[size_position(1),size_position(2),min(size(logo,2),round(size_position(3)*.8)),min(size(logo,2),round(size_position(3)*.8))/size(logo,2)*size(logo,1)]);
    set(s,'Position',[0,0,1,1])
    overview_fontsize = (min(size(logo,2),round(size_position(3)*.8))/size(logo,2)*size(logo,1))/30;
    set(overview,'visible','on')
    pause(1)
    %---------------------------------------------------------------------%
    
    
    
    %% define directory
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.35,.5,.05],'String','loading parameters...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start_total = now();
    %---------------------------------------------------------------------%
    
    
    if read_from_xlsx ==1
        % determine directory from txt- and excel-file
        name_file = info_input{2};
        cd(name_directory)
        
        % read excel-files
        xlsx_name = 'copy_parameters.xlsx';
        cd(name_directory)
        copyfile(name_file,xlsx_name)
        cd(name_directory)
        tbl_input_data = readtable(xlsx_name,'sheet','required_input_data');
        tbl_parameters = readtable(xlsx_name,'sheet','parameters');
        tbl_raster_selected = readtable(xlsx_name,'sheet','raster_selected');
        delete(xlsx_name)
        
        path_result = char(tbl_input_data{2,2});
        raster_selected.name = char(tbl_raster_selected.name(1));
        raster_selected.legend = tbl_raster_selected.legend';
        raster_selected.color = [tbl_raster_selected.color_r,tbl_raster_selected.color_g,tbl_raster_selected.color_b];
        
    else
        % alternative: define path manually
        path_result = 'DIRECTORY IN WHICH THE RESULT IS SAVED';
        raster_selected.name = 'landuse';
        raster_selected.legend = {'pasture','cropland','forest','track','road','sealed'};
        raster_selected.color = [[86, 129, 45]/255;
            [144, 131, 12]/255;
            [43, 64, 27]/255;
            [102, 105, 99]/255;
            [177, 108, 37]/255;
            [141, 36, 37]/255];
    end
    
    %% define requested output
    
    % definition, which output is saved
    save_grids = 0;
    save_river_points = 0;
    save_dam_points = 0;
    save_basins_selected = 0;
    save_basins_as_ascii = tbl_parameters{25,2};
    save_curves_as_excel = tbl_parameters{26,2};
    
    save_options = struct('save_grids',save_grids,...
        'save_river_points',save_river_points,...
        'save_dam_points',save_dam_points,...
        'save_basins_selected',save_basins_selected,...
        'save_basins_as_ascii',save_basins_as_ascii,...
        'save_curves_as_excel',save_curves_as_excel);
    
    
    
    %% define requested plots
    
    % definition, which output is plotted
    plot_exitcodes = tbl_parameters{27,2};
    plot_spatial_overview = tbl_parameters{28,2};
    plot_factsheet_p1 = tbl_parameters{29,2};
    plot_factsheet_p2 = tbl_parameters{30,2};
    plot_dam_comparison = tbl_parameters{31,2};
    plot_curve_comparison = tbl_parameters{32,2};
    plot_visibility = tbl_parameters{33,2};
    
    plot_options = struct('plot_exitcodes',plot_exitcodes,...
        'plot_spatial_overview',plot_spatial_overview,...
        'plot_factsheet_p1',plot_factsheet_p1,...
        'plot_factsheet_p2',plot_factsheet_p2,...
        'plot_dam_comparison',plot_dam_comparison,...
        'plot_curve_comparison',plot_curve_comparison,...
        'plot_visibility',plot_visibility);
    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start_total)*24*60;
    fprintf('Parameter setting completed: %3.3f min\n',duration);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.35,.5,.05],'String',sprintf('finished: %3.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %---------------------------------------------------------------------%
    
    
    %% load data
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.3,.5,.05],'String','loading data...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start = now();
    %---------------------------------------------------------------------%
    
    
    cd(path_result)
    load('input_grids_used.mat')
    load('river_points.mat')
    load('basins_selected.mat')
    
    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start)*24*60;
    fprintf('Data import completed: %3.3f min\n',duration);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.3,.5,.05],'String',sprintf('finished: %3.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %---------------------------------------------------------------------%
 
    
    %% export requested results
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.25,.5,.05],'String','export results...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start = now();
    %---------------------------------------------------------------------%
    
    
    %%% save results %%%
    fun_9_summary_and_save_output(save_options)
    
    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start)*24*60;
    fprintf('Output saved: %3.3f min\n',duration);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.25,.5,.05],'String',sprintf('finished: %3.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
    end
    %---------------------------------------------------------------------%
    
    
    %% plot requested figures
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.1,.2,.5,.05],'String','plotting results...','linestyle','none','fontsize',overview_fontsize)
        pause(1)
    end
    %-- calculate duration -----------------------------------------------%
    time_start = now();
    %---------------------------------------------------------------------%
    
    
    %%% plot figures %%%
    fun_10_plot_results(plot_options,raster_selected);

    
    %-- calculate duration -----------------------------------------------%
    time_end = now;
    duration = (time_end-time_start)*24*60;
    duration_total = (time_end-time_start_total)*24*60;
    fprintf('Figures saved: %3.3f min\n',duration);
    fprintf('Total duration: %3.3f min\n',duration_total);
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.4,.2,.5,.05],'String',sprintf('finished: %3.3f min',duration),'linestyle','none','fontsize',overview_fontsize)
        pause(1)
        annotation('textbox',[.45,.12,.5,.05],'String',sprintf('Total duration: %3.3f min\n',duration_total),'linestyle','none','fontsize',overview_fontsize)
    end
    %---------------------------------------------------------------------%
    
    diary off
    
catch ME
    
    %-- update status figure ---------------------------------------------%
    if ishandle(overview)==1
        set(0,'CurrentFigure',overview)
        annotation('textbox',[.45,.12,.4,.05],'String','ERROR. Read logfile for more details.','linestyle','none','fontsize',overview_fontsize)
    end
    %---------------------------------------------------------------------%
    
    rethrow(ME)
    diary off
    
end
