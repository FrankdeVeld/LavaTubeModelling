%% Goal of this function: export a station file grid for IGMAS
% Author: Frank de Veld
% Date: June 17th, 2020
% Output:  one .station file

% Function to also export a stations file. There is no elevation assumed,
% and also just a standard grid. While it is interesting to vary station
% shapes, it is advised to just define a fine grid here and model with
% this, and afterwards do operations on the obtained data (e.g. different
% patterns)

% This function is called automatically when running the Voxel_main.m
% script
function Export_stations_file(Model_limits, Station_spacing, Model_name)
    feature('DefaultCharacterSet', 'UTF8'); 
    File_name = [Model_name,'.stations'];
    fid = fopen(File_name,'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n'); % Result in proper encoding
    Temp_text = ['<geodata> \n' ...
    '<projection name="unknown" units="m"></projection> \n']; % Open the geodata class
    fprintf(fid, Temp_text);
    x_increment = Station_spacing;
    y_increment = Station_spacing;
    for a=1:(round((Model_limits(2)-Model_limits(1))/x_increment)+1)
        for b=1:(round((Model_limits(4)-Model_limits(3))/y_increment)+1)
            Temp_text = ['<vertex x="%4.2f" y="%4.2f" z="%4.2f"> \n</vertex> \n'];
            Temp_text = sprintf(Temp_text,Model_limits(1) + x_increment*(a-1),Model_limits(3) + y_increment*(b-1),Model_limits(6));
            fprintf(fid,Temp_text);
        end
    end
    clear a b
    fprintf(fid, '</geodata>');
end