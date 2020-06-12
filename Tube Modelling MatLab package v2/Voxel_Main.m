%% Script for making density-based cavity models for IGMAS modelling
% Authors: Frank de Veld
% Date: June 10th, 2020
% Output: .model file, .vxo file x2, .stations file
% Version: 2.0

% v 1.0:    May 25th, 2020. First, working version
% v 2.0:    June 9th, 2020. major cut in running time and big improvement in
%           efficiency. Change to 'inShape' function for checking whather points
%           are in a cavity. Script has been made clearer and smaller
% v 2.1:    June 11th, 2020. Made cavity density variable, added a way of
%           visualisation. Added importing of basic shapes. Todo: add error
%           messages, make importing more variable

% Aim of this script: Go from a coordinate-based file to an IGMAS modelling
% file

% Advantages: intuitive to model, all possible caves can be modelled,
% script is very understandable 
% Disadvantage: voxels must be placed strictly on a grid (so keep enough distance between voxels)

% Running the script should take about 0-20 seconds (estimated)

% What can be changed?
% - Reducing factor: how many data points do you find acceptable (Step 0)
% - Sizing factor: how much larger you want your model to be, compared
% to its original size (Step 0)
% - Horizontal stretching factor: how much do you want your cavity to be
% stretched horizontally (Step 0)

% - Cavity depth: The distance between cavity and stations (Step 1)
% - Rock density: what the density of the surrounding rock is (Step 1)
% - Cavity density: density of interior of cavity (so, probably 0) (Step 1)
% 
% - Name of the model (Step 1)
% - Edge factor: how much extra space you want in the x- and y-directions
% to prevent boundary effects (Step 1)

% - Station resolution: how many stations do you want on the grid (Step 3)

clear all
close all
%% Step 0: Import the coordinates and reduce the resolution

% The file either is an csv-file (use importdata) or a .stl file (use stlread)

% 'Cave_info' consists of faces and vertices
% 'stlread' is a function from the MatLab file exchange, see https://nl.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader
Cave_info = stlread('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Lava tube models\galapagos-lava-tube\source\Mesh SCruz2 levigata\Mesh SCruz2 levigata.stl');
Cavity_coordinates = Cave_info.vertices;

% Other option: make a standard shape. Options: 'Cylinder', 'Sphere' and
% 'Beam'
% % 
% Cylinder_length = 100; % Horizontal length in meter
% Cylinder_radius = 15; % Vertical radius in meter
% Cavity_coordinates = Import_shape('Cylinder_H', Cylinder_length, 0, 0, Cylinder_radius);

Cylinder_height = 100; % Horizontal length in meter
Cylinder_radius = 15; % Vertical radius in meter
Cavity_coordinates = Import_shape('Cylinder_V', 0, 0, Cylinder_height, Cylinder_radius);

% Sphere_radius = 10;
% Cavity_coordinates = Import_shape('Sphere', 0, 0, 0, Sphere_radius);

% Beam_length = 60;
% Beam_width = 30;
% Beam_heigth = 10;
% Cavity_coordinates = Import_shape('Beam', Beam_length, Beam_width, Beam_heigth, 0);

% If the number of points is much more than you need, you can reduce the
% number with this parameter:
Reduce_factor = 1;

% It takes some time and right now we don't use it, so it is commented out
% Cavity_coordinates = Reduce_coordinates(Cavity_coordinates, Reduce_factor);

% What this does; make the cavity Sizing_factor times larger, or stretch it
% by a factor Stretching_factor_h horizontally.
% Vertical stretching is a bit trickier, or shape-preserving stretching
% Alternative; stretch it until the average radius/length is .../... meter
Sizing_factor = 1;
Stretching_factor_h = 1;

Cavity_coordinates = Sizing_factor*Cavity_coordinates;
Cavity_coordinates = Stretching_factor_h*Cavity_coordinates;

%% Step 1: Define parameters and setup 

% A grid will be made, and here you can choose the spacing between points
% (same in all directions)
% It is multiplied with the sizing factor for the cavity automatically
Grid_spacing = 5;
Grid_spacing = Grid_spacing*Sizing_factor;
Rock_density = 0; % density in g/cm³
Cavity_density = 2.5; % density in g/cm³
Model_name = ['Cylinder_V_gs5_verification_test']; % The name you want your model to have
Edging_factor = 0.50; % How much additional rock you want around the cavity, in order to prevent edge effects to play a dominant role. 
% This can be set 0 without too many problems, but from about 0.25 on you
% probably won't see the edge effects anymore
Cavity_depth = 10; % How much space filled with rock there is between the top of the cavity and the stations (on the ground)

Gridded_coords = {};
Sizing_addition = [];
Model_limits = [];
Cavity_limits = [];
for a=1:3 % In this loop (3 or 3x2 coordinates) the grid is created and the limits for the cavity and the model (with edging factor) are calculated
    Gridded_coords(a) = {round(min(Cavity_coordinates(:,a))):Grid_spacing:(round(max(Cavity_coordinates(:,a)))-Grid_spacing)};
    Sizing_addition(a) = round((Gridded_coords{1,a}(end)-Gridded_coords{1,a}(1))*Edging_factor/Grid_spacing)*Grid_spacing;
    Cavity_limits(2*a-1) = round(min(Cavity_coordinates(:,a)));
    Cavity_limits(2*a) = round(max(Cavity_coordinates(:,a)));
    if a<3
        Model_limits(2*a-1) = round(min(Cavity_coordinates(:,a)))-Sizing_addition(a);
        Model_limits(2*a) = round(max(Cavity_coordinates(:,a)))+Sizing_addition(a);
    else
        Model_limits(2*a-1) = round(min(Cavity_coordinates(:,a)))-Cavity_depth;
        Model_limits(2*a) = round(max(Cavity_coordinates(:,a)))+Cavity_depth;
    end
end

for b=1:3
    Gridded_coords(b) = {Model_limits(2*b-1):Grid_spacing:(Model_limits(2*b)-Grid_spacing)};
end % Update the grid to include the whole model, including edge effect sizing

% Shift everything lower/higher with the maximum z-coordinate of the model,
% such that the model ends at z=0, the cavity begins at z=0-Cavity_depth
% and all limits are changed accordingly
Gridded_coords{1,3} = Gridded_coords{1,3} - Model_limits(6);
Cavity_coordinates(:,3) = Cavity_coordinates(:,3) - Model_limits(6);
Cavity_limits(5) = Cavity_limits(5) - Model_limits(6);
Cavity_limits(6) = Cavity_limits(6) - Model_limits(6);
Model_limits(5) = Model_limits(5) - Model_limits(6);
Model_limits(6) = Model_limits(6) - Model_limits(6);

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
Visualisation_bool = true; % Also export a file with 'just the cavity'
File_name = [Model_name,'_With_Cavity.vxo']; % The model name defined earlier is added
Cavity_bool = true; % This is to use the following function twice, instead of making separate functions for it
Voxel_Creation(Rock_density,Cavity_density,Cavity_coordinates,Grid_matrix_3n,Cavity_bool,Model_limits,Grid_spacing,File_name,Visualisation_bool)

% Without cavity
Visualisation_bool = false;
File_name = [Model_name,'_No_Cavity.vxo'];
Cavity_bool = false;
Voxel_Creation(Rock_density,Cavity_density,Cavity_coordinates,Grid_matrix_3n,Cavity_bool,Model_limits,Grid_spacing,File_name,Visualisation_bool)


% Export model file
Write_basic_model_file(Rock_density, Cavity_density, Model_name, Model_limits, Cavity_limits, length(Gridded_coords{1,2})-1) 
%% Step 3: Define the stations
Station_resolution = 50; % Advised is to have the number of stations per m² be comparable to the number of voxels per m³, or at least in the same order of magnitude

Export_stations_file(Model_limits,Station_resolution,Model_name) % ASSUMING flat land, also currently only a set grid is created without other patterns

%% Step 4: Import the model file, the voxel files and the station file in IGMAS, do the anomaly calculation and export the .stations files

% In order to make Step 5 work well, only calculate 'gz' and not the other
% signals. This is the interesting signal

%% Step 5: Compare the calculated anomalies

% There is a different file called 'Data_handling.m' for this right now.
% Since the whole modelling process can't be automated right now, these
% scripts are better separated.
