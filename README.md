# LOCASIN

The open-source tool *LOCASIN* enables an automated and rapid detection, characterization and evaluation of basin locations. 
The coding was done in MATLAB R2018b. 
The analysis includes the determination of the optimal dam axis orientation, the dam geometry, the basin area and the basin volume. 

## (1) Program Structure

[LOCASIN.m](./code/LOCASIN.m)

1.  Preprocessing
    *  [define_input_directory_and_file.txt](./code/define_input_directory_and_file.txt)
    *  [USER_INPUT.m](./code/USER_INPUT.m) / [USER_INPUT_from_xlsx.m](./code/USER_INPUT_from_xlsx.m)
    *  [fun_1_river_analysis.m](./code/fun_1_river_analysis.m)
        *  [fun_1n1_ASCIimport.m](./code/fun_1n1_ASCIimport.m)
2.  Analysis of the Dam Positioning
    *  [fun_2_determine_shortest_dam.m](./code/fun_2_determine_shortest_dam.m)
3.  Basin Analysis
    *  [fun_3_determine_river_points_in_basin.m](./code/fun_3_determine_river_points_in_basin.m)
    *  [fun_4_determine_dam_characteristics.m](./cdoe/fun_4_determine_dam_characteristics.m)
    *  [fun_5_determine_basin_area_and_volume.m](./code/fun_5_determine_basin_area_and_volume.m)
    *  [fun_6_evaluate_basins.m](./code/fun_6_evaluate_basins.m)
4.  Basin Combination
    *  [fun_7_selection_of_basin_combination.m](./code/fun_7_selection_of_basin_combination.m)
    *  [fun_8_determine_depth_storage_area_curves.m](./code/fun_8_determine_depth_storage_area_curves.m)
5.  Postprocessing
    *  [fun_9_summary_and_save_output.m](./code/fun_9_summary_and_save_output.m)
        *  [fun_9n1_save_basins_as_ascii.m](./code/fun_9n1_save_basins_as_ascii.m)
        *  [fun_9n2_save_basins_characteristics_as_xlsx.m](./code/fun_9n2_save_basins_characteristics_as_xlsx.m)
    *  [fun_10_plot_results.m](./code/fun_10_plot_results.m)
        *  [plot_1_spatial_exit_codes.m](./code/plot_1_spatial_exit_codes.m)
        *  [plot_2_spatial_all_basins.m](./code/plot_2_spatial_all_basins.m)
        *  [plot_3_single_basin_characteristics.m](./code/plot_3_single_basin_characteristics.m)
            *  [plot_3n1_spatial_basin.m](./code/plot_3n1_spatial_basin.m)
            *  [plot_3n2_updateContours.m](./code/plot_3n2_updateContours.m)
        *  [plot_4_potential_dam_characteristics_portrait.m](./code/plot_4_potential_dam_characteristics_portrait.m)
        *  [plot_5_curves_depth_storage_area.m](./code/plot_5_curves_depth_storage_area.m)
   

<br>


## (2) Additional Files

1.  [LOCASIN_plotting.m](./code/LOCASIN_plotting.m)
    *  [define_input_directory_and_file_plotting.txt](./code/define_input_directory_and_file_plotting.txt)
2.  [fun_replace1string.m](./code/fun_replace1string.m)
3.  Figure
    *  [locasin.jpg](./code/locasin.jpg)
    *  [locasin_plotting.jpg](./code/locasin_plotting.jpg)
    *  [locasin_icon.jpg](./code/locasin_icon.jpg)


<br>

## (3) Standalone Application

1. [Installer](./standalone_applocation/installer)
    *  Installer_mrc.exe
    *  Installer_web.exe
2. LOCASIN
3. LOCASINplotting


<br>

## (4) Testcase

