%% Goal of this function: provide the outline of the cavity to the applet
% Author: Frank de Veld
% Date: June 22th, 2020
% Output: A matrix of outline coordinates

Lava_tube_model_array = {'Name1','Name2','Name3','Name4'};
for a=1:4
    Lava_tube = Lava_tube_model_array{a};

    if (strcmp(Lava_tube,'Name1')==true)
        File_directory = '.\Model files\Name1.stl';
    elseif (strcmp(Lava_tube,'Name2')==true)
        File_directory = '.\Model files\Name2.stl';
    elseif (strcmp(Lava_tube,'Name3')==true) 
        File_directory = '.\Model files\Name3.stl';
    elseif (strcmp(Lava_tube,'Name4')==true)
        File_directory = '.\Model files\Name4.stl';
    end          
    
    Cave_info = stlread(File_directory);
    Coordinates = Cave_info.vertices;
    Outline = boundary(Coordinates(:,1),Coordinates(:,2));
    
    x_coords = Coordinates(:,1);
    y_coords = Coordinates(:,2);
    Outline_coordinates = [x_coords(Outline),y_coords(Outline)];
    Outline_name = ['Outline_',Lava_tube];
    save(Outline_name,'Outline_coordinates');
end
