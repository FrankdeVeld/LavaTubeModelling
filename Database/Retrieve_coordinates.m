%% Goal of this function: provide the coordinates of the cavity to the applet
% Author: Frank de Veld
% Date: June 17th, 2020
% Output: A matrix of Cavity coordinates 

% A bit of an overkill, as we just would like the outline really

% This function is called automativally when the 'Show cavity' box is
% checked in the app 'Data_visualisation_app'

function [Coordinates] = Retrieve_coordinates(Lava_tube)

    if (strcmp(Lava_tube,'LAGOS')==true)
        File_directory = 'C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Lava tube models\cueva-de-los-siete-lagos\source\agua_text\Los lagos ok.stl';
    elseif (strcmp(Lava_tube,'GALA')==true)
        File_directory = 'C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Lava tube models\galapagos-lava-tube\source\Mesh SCruz2 levigata\Mesh SCruz2 levigata.stl';
    elseif (strcmp(Lava_tube,'JAMEO')==true) 
        File_directory = 'C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Lava tube models\lava-tube-jameo-de-la-prendes\source\1 Jameo de la Prendes\1 Jameo de la Prendes.stl';
    elseif (strcmp(Lava_tube,'ETNA')==true)
        File_directory = 'C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Lava tube models\lavatube-mount-etna\source\2ae044b5e89b469a9ac67baea1627996.fbx';
    end       
    Cave_info = stlread(File_directory);
    Coordinates = Cave_info.vertices;
end
