# LOCASIN

The open-source tool *LOCASIN* enables an automated and rapid detection, characterization and evaluation of basin locations. 
The coding was done in MATLAB R2018b. 
The analysis includes the determination of the optimal dam axis orientation, the dam geometry, the basin area and the basin volume. 

## (1) Program Structure

[LOCASIN.m](./code/LOCASIN.m)<br>
LOCation detection of retention and detention bASINs

1.  Preprocessing
    *  [define_input_directory_and_file.txt](./code/define_input_directory_and_file.txt):<br>
         definition of the directory where the excel-user-input is located and name of the excel-file
    *  [USER_INPUT.m](./code/USER_INPUT.m) / [USER_INPUT_from_xlsx.m](./code/USER_INPUT_from_xlsx.m):<br>
         summary of all user inputs which are required for the determination of basin locations
    *  [fun_1_river_analysis.m](./code/fun_1_river_analysis.m):<br>
         function to analyze and characterize river points
        *  [fun_1n1_ASCIimport.m](./code/fun_1n1_ASCIimport.m):<br>
            import ASCII raster data into matlab including the coordinates of the reference system and the size information
2.  Analysis of the Dam Positioning
    *  [fun_2_determine_shortest_dam.m](./code/fun_2_determine_shortest_dam.m):<br>
         function to determine the dam orientation with the shortest dam length
3.  Basin Analysis
    *  [fun_3_determine_river_points_in_basin.m](./code/fun_3_determine_river_points_in_basin.m):<br>
         function to determine all potential points in the basin to estimate the reqired extent of the raster for the basin calculation
    *  [fun_4_determine_dam_characteristics.m](./cdoe/fun_4_determine_dam_characteristics.m):<br>
         function to determine the geometry of the dam (volume, height, ...)
    *  [fun_5_determine_basin_area_and_volume.m](./code/fun_5_determine_basin_area_and_volume.m):<br>
         function to determine the inundation area of the basin and the subsequent depth distribution and volume
    *  [fun_6_evaluate_basins.m](./code/fun_6_evaluate_basins.m):<br>
         function to evaluate and exclude dam heights based on the dam and basin characteristics
4.  Basin Combination
    *  [fun_7_selection_of_basin_combination.m](./code/fun_7_selection_of_basin_combination.m):<br>
         function to select the best basin combination based on the evaluation cirteria and overlapping flooding areas
    *  [fun_8_determine_depth_storage_area_curves.m](./code/fun_8_determine_depth_storage_area_curves.m):<br>
         function to calculate curves of the relation between water depth, storage volume and flooding area for the basins of the selected combination
5.  Postprocessing
    *  [fun_9_summary_and_save_output.m](./code/fun_9_summary_and_save_output.m):<br>
         function to fummary the output and to save the results based on the user definitions (selected detail of the output)
        *  [fun_9n1_save_basins_as_ascii.m](./code/fun_9n1_save_basins_as_ascii.m):<br>
            function to save the basin characteristics as ascii raster file which can be imported to GIS programs (dam heigths are defined positive, while water depths are defined negative)
        *  [fun_9n2_save_basins_characteristics_as_xlsx.m](./code/fun_9n2_save_basins_characteristics_as_xlsx.m):<br>
            function to save curves (depth-storage-area) of the basins in an excel-file
    *  [fun_10_plot_results.m](./code/fun_10_plot_results.m):<br>
         function to plot the results based on predefined standard plots for the visual representation of the most important results
        *  [plot_1_spatial_exit_codes.m](./code/plot_1_spatial_exit_codes.m):<br>
            plot to analyze the procedure of the river analysis and basin selection
        *  [plot_2_spatial_all_basins.m](./code/plot_2_spatial_all_basins.m):<br>
            plot show the locations of the dams and basins of the selected basin combination including dam heights and water depths
        *  [plot_3_single_basin_characteristics.m](./code/plot_3_single_basin_characteristics.m):<br>
            fact sheet, page 1: plot of the dam and basin characteristics for the best dam height
            *  [plot_3n1_spatial_basin.m](./code/plot_3n1_spatial_basin.m):<br>
               function to plot one or more dams and basins in a top view including background maßs and elevation lines
            *  [plot_3n2_updateContours.m](./code/plot_3n2_updateContours.m)
        *  [plot_4_potential_dam_characteristics_portrait.m](./code/plot_4_potential_dam_characteristics_portrait.m):<br>
            fact sheet, page 2: plot the dam characteristics, basin characteristics and evaluation criteria for all potential dam heights for one or more dam sites
        *  [plot_5_curves_depth_storage_area.m](./code/plot_5_curves_depth_storage_area.m):<br>
            plot curves of water depths, storage volume and flooding area to compare multiple dam sites
   

<br>


## (2) Additional Files

1.  [LOCASIN_plotting.m](./code/LOCASIN_plotting.m):<br>
      separate plotting of the figures form already calculated *.mat-results
    *  [define_input_directory_and_file_plotting.txt](./code/define_input_directory_and_file_plotting.txt):<br>
         definition of the directory where the excel-user-input is located and name of the excel-file
2.  [fun_replace1string.m](./code/fun_replace1string.m):<br>
      function to replace two strings in a txt-file, e.g. comma by point
3.  Figure
    *  [locasin.jpg](./code/locasin.jpg):<br>
         splash screen of the LOCASIN tool
    *  [locasin_plotting.jpg](./code/locasin_plotting.jpg):<br>
         splash screen of the LOCASINplotting tool
    *  [locasin_icon.jpg](./code/locasin_icon.jpg):<br>
         icon of the LOCASIN tool

<br>

## (3) Standalone Application

1. [Installer](./standalone_application/installer)
    *  Installer_web.exe: <br>
         installer for the MATLAB version for standalone applications (data is downloaded form the web)
2. [LOCASIN](./standalone_application/LOCASIN)
    *  LOCASIN.exe:<br>
         standalone application of LOCASIN
    *  [define_input_directory_and_file.txt](./code/define_input_directory_and_file.txt):<br>
         definition of the directory where the excel-user-input is located and name of the excel-file
3. [LOCASINplotting](./standalone_application/LOCASINplotting)
    *  LOCASINplotting.exe:<br>
         standalone application of LOCASINplotting
    *  [define_input_directory_and_file.txt](./code/define_input_directory_and_file_plotting.txt):<br>
         definition of the directory where the excel-user-input is located and name of the excel-file

<br>

## (4) Testcase

1. User Information

2. Input Data

3. Results

