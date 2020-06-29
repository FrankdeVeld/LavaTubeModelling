%% Goal of this script: making density-based cavity models for IGMAS modelling
% Author: Frank de Veld
% Date: June 18th, 2020
% Output: .model file, .vxo file x2, .stations file
% Version: 2.0

% v 1.0:    May 25th, 2020. First, working version
% v 2.0:    June 9th, 2020. Major cut in running time and big improvement in
%           efficiency. Change to 'inShape' function for checking whather points
%           are in a cavity. Script has been made clearer and smaller
% v 2.1:    June 11th, 2020. Made cavity density variable, added a way of
%           visualisation. Added importing of basic shapes.  Made file importing more variable 
%           Sdded error messages 
% v 2.5:    June 15th, 2020. No export anymore of 'cavity_less' file. No
%           voxels created outside of cavity. Cavity has density of density
%           difference. Stations account for edge spacing. no visualisation
%           file

% Aim of this script: Go from a coordinate-based file to an IGMAS modelling
% file

% Advantages: intuitive to model, all possible caves can be modelled,
% script is very understandable 
% Disadvantage: voxels must be placed strictly on a grid (so keep enough distance between voxels)

% Running the script should take about 0-20 seconds (estimated)

% - Station resolution: how many stations do you want on the grid (Step 3)

% clear all
% close all

% These are the parameters that need to be changed constantly for making
% the database
Lava_tube_model_array = {'Name1','Name2','Name3','Name4'};
Model_multiplication_array = [1, 2, 5, 10]; % Leads to cavities of about 20 m diameter to 200 m diameter
Model_density_array = [1.5, 2.25, 3]; % Densities are easy to interpret and interpolate
Cavity_depth_array = [10, 20, 40, 80]; % Such that most cavities are still structurally stable

Lava_tube = Lava_tube_model_array{2};

Create_voxels(Lava_tube,1,1.5,10,2.5,2.5)

for a=3:4
    for b=1:3
        for c=1:4
            Model_multiplication = Model_multiplication_array(a);
            Model_density = Model_density_array(b);
            Cavity_depth = Cavity_depth_array(c);
            Station_spacing = 2.5*Model_multiplication; % Advised is to have the number of stations per m² be comparable to the number of voxels per m³, or at least in the same order of magnitude
            Grid_spacing = 2.5*Model_multiplication; % Advised is to have the number of stations per m² be comparable to the number of voxels per m³, or at least in the same order of magnitude

            Create_voxels(Lava_tube,Model_multiplication,Model_density,Cavity_depth,Station_spacing,Grid_spacing)
        end
    end
end

function Create_voxels(Lava_tube,Model_multiplication,Model_density,Cavity_depth,Station_spacing,Grid_spacing)

    if (strcmp(Lava_tube,'Name1')==true)
        File_directory = '.\Model files\Name1.stl';
    elseif (strcmp(Lava_tube,'Name2')==true)
        File_directory = '.\Model files\Name2.stl';
    elseif (strcmp(Lava_tube,'Name3')==true) 
        File_directory = '.\Model files\Name3.stl';
    elseif (strcmp(Lava_tube,'Name4')==true)
        File_directory = '.\Model files\Name4.stl';
    end     

    Model_name = strcat(Lava_tube,'_MM_',num2str(Model_multiplication),'_MD_',num2str(Model_density),'_CD_',num2str(Cavity_depth));
    %% Step 0: Import the coordinates and reduce the resolution

    % 'Cave_info' consists of faces and vertices
    % 'stlread' is a function from the MatLab file exchange, see https://nl.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader
    Cave_info = stlread(File_directory);
    Cavity_coordinates = Cave_info.vertices;

    % What this does; make the cavity Sizing_factor times larger to make them
    % more fit for the lunar environment
    Sizing_factor = Model_multiplication;

    Cavity_coordinates = Sizing_factor*Cavity_coordinates;

    %% Step 1: Define parameters and setup 

    % A grid will be made, and here you can choose the spacing between points
    % (same in all directions)
    % It is multiplied with the sizing factor for the cavity automatically
    % Grid_spacing = 5;
    % Grid_spacing = Grid_spacing*Sizing_factor;
    % Rock_density = 2.5; % density in g/cm³
    % Cavity_density = 0; % density in g/cm³
    % Model_density = abs(Rock_density - Cavity_density);

    Model_name = [Model_name]; % The name you want your model to have
    Edging_factor = 0.50; % How much additional spacing you want around the cavity, in order to prevent edge effects to play a dominant role. 
    % This can be set 0 without too many problems, but from about 0.25 on you
    % probably won't see the edge effects anymore. This does not affect the
    % performance anymore, and just places extra stations
    % Cavity_depth = 10; % How much space filled with rock there is between the top of the cavity and the stations (on the ground)

    if Grid_spacing <=0
        error('Error: grid spacing must be positive')
    end
    if ( (max(Cavity_coordinates(:,1))-min(Cavity_coordinates(:,1)))/Grid_spacing * (max(Cavity_coordinates(:,2))-min(Cavity_coordinates(:,2)))/Grid_spacing * (max(Cavity_coordinates(:,3))-min(Cavity_coordinates(:,3)))/Grid_spacing > 10000000) 
        warning('Warning: many voxels present (> 10 000 000), modelling process can take a while')
    elseif  ( (max(Cavity_coordinates(:,1))-min(Cavity_coordinates(:,1)))/Grid_spacing * (max(Cavity_coordinates(:,2))-min(Cavity_coordinates(:,2)))/Grid_spacing * (max(Cavity_coordinates(:,3))-min(Cavity_coordinates(:,3)))/Grid_spacing > 30000000) 
        error('Error: too many voxels present (> 30 000 000). Reduce by adjusting the grid spacing')
    end

    if( (max(Cavity_coordinates(:,1))-min(Cavity_coordinates(:,1)))/Grid_spacing < 1 || (max(Cavity_coordinates(:,2))-min(Cavity_coordinates(:,2)))/Grid_spacing < 1 || (max(Cavity_coordinates(:,3))-min(Cavity_coordinates(:,3)))/Grid_spacing < 1)
        error('Error: grid spacing too high; zero voxels detected in one of the directions')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Gridded_coords = {};
    Sizing_addition = [];
    Model_limits = [];
    Cavity_limits = [];
    for a=1:3 % In this loop (3 or 3x2 coordinates) the grid is created and the limits for the cavity and the model (with edging factor) are calculated
        Gridded_coords(a) = {round(min(Cavity_coordinates(:,a))):Grid_spacing:(round(max(Cavity_coordinates(:,a)))-Grid_spacing)};
        Sizing_addition(a) = round((Gridded_coords{1,a}(end)-Gridded_coords{1,a}(1))*Edging_factor/Grid_spacing)*Grid_spacing;
        Cavity_limits(2*a-1) = round(min(Cavity_coordinates(:,a)));
        Cavity_limits(2*a) = round(max(Cavity_coordinates(:,a)));
        if a==1
            Model_limits(2*a-1) = round(min(Cavity_coordinates(:,a)))-Sizing_addition(a);
            Model_limits(2*a) = round(max(Cavity_coordinates(:,a)))+Sizing_addition(a);
        elseif a ==2
            Model_limits(2*a-1) = round(min(Cavity_coordinates(:,a)))-Sizing_addition(a); % Can be + for 'cut-off' cavity
            Model_limits(2*a) = round(max(Cavity_coordinates(:,a)))+Sizing_addition(a); % Can be - for 'cut-off' cavity
        elseif a==3
            Model_limits(2*a-1) = round(min(Cavity_coordinates(:,a)))-Cavity_depth;
            Model_limits(2*a) = round(max(Cavity_coordinates(:,a)))+Cavity_depth;
        end
    end

    % for b=1:3
    %     Gridded_coords(b) = {Model_limits(2*b-1):Grid_spacing:(Model_limits(2*b)-Grid_spacing)};
    % end % Update the grid to include the whole model, including edge effect sizing

    % Shift everything lower/higher with the maximum z-coordinate of the model,
    % such that the model ends at z=0, the cavity begins at z=0-Cavity_depth
    % and all limits are changed accordingly
    Gridded_coords{1,3} = Gridded_coords{1,3} - Model_limits(6);
    Cavity_coordinates(:,3) = Cavity_coordinates(:,3) - Model_limits(6);
    Cavity_limits(5) = Cavity_limits(5) - Model_limits(6);
    Cavity_limits(6) = Cavity_limits(6) - Model_limits(6);
    Model_limits(5) = Model_limits(5) - Model_limits(6);
    Model_limits(6) = Model_limits(6) - Model_limits(6);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This is still to create a grid, but in the format of an n x n x n
    % meshgrid. This makes it easier to create 'logicals', which speed uo the
    % process later on

    % Not the addition of half a grid spacing! This is since voxle coordinates
    % actually signifiy the locations of the centers of voxels.
    [X,Y,Z] = meshgrid(Gridded_coords{1,1}+Grid_spacing/2,Gridded_coords{1,2}+Grid_spacing/2,Gridded_coords{1,3}+Grid_spacing/2);
    Logical = true(length(Gridded_coords{1,2}),length(Gridded_coords{1,1}),length(Gridded_coords{1,3}));
    % The Grid_matrix_3n is an nx3 matrix consisting of every coordinate
    % considered in the model
    Grid_matrix_3n = [X(Logical),Y(Logical),Z(Logical)];

    %% Step 2: Export a voxel file with cavity, one without and a basic .model file
    % Now that the grid has been created and is in the right format and the limits
    % are known, the voxel models can be made

    % 'Voxel_Creation' is a separate script written by the author and can be
    % found in the same file set
    % With cavity
    File_directory = [Model_name,'.vxo']; % The model name defined earlier is added
    Cavity_bool = true; % This is to use the following function twice, instead of making separate functions for it
    Voxel_creation(Model_density,Cavity_coordinates,Grid_matrix_3n,File_directory)

    % Export model file
    Write_basic_model_file(Model_density, Model_name, Model_limits, Cavity_limits, length(Gridded_coords{1,2})-1) 
    %% Step 3: Define the stations

    Export_stations_file(Model_limits,Station_spacing,Model_name) % ASSUMING flat land, also currently only a set grid is created without other patterns

end
%% Step 4: Import the model file, the voxel files and the station file in IGMAS, do the anomaly calculation and export the .stations files

% In order to make Step 5 work well, only calculate 'gz' and not the other
% signals. This is the interesting signal in the first place

%% Step 5: Compare the calculated anomalies

% There is a different file called 'Data_handling.m' for this right now.
% Since the whole modelling process can't be automated right now, these
% scripts are better separated.
