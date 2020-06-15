% Function to also export a stations file. There is no elevation assumed,
% and also just a standard grid. While it is interesting to vary station
% shapes, it is advised to just define a fine grid here and model with
% this, and afterwards do operations on the obtained data (e.g. different
% patterns)
function Export_stations_file(Station_limit_vector,Station_resolution, Model_name)
    feature('DefaultCharacterSet', 'UTF8'); 
    File_name = [Model_name,'_Stations.stations'];
    fid = fopen(File_name,'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n'); % Result in proper encoding
    Temp_text = ['<geodata> \n' ...
    '<projection name="unknown" units="m"></projection> \n']; % Open the geodata class
    fprintf(fid, Temp_text);
    x_increment = (Station_limit_vector(2)-Station_limit_vector(1))/Station_resolution;
    y_increment = (Station_limit_vector(4)-Station_limit_vector(3))/Station_resolution;
    for a=1:(Station_resolution+1)
        for b=1:(Station_resolution+1)
            Temp_text = ['<vertex x="%4.2f" y="%4.2f" z="%4.2f"> \n</vertex> \n'];
            Temp_text = sprintf(Temp_text,Station_limit_vector(1) + x_increment*(a-1),Station_limit_vector(3) + y_increment*(b-1),Station_limit_vector(6));
            fprintf(fid,Temp_text);
        end
    end
    clear a b
    fprintf(fid, '</geodata>');
end