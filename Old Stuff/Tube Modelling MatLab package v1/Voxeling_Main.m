%% Script for making density-based cavity models for IGMAS modelling
% Authors: Frank de Veld
% Date: May 25th, 2020
% Output: .model file, .vxo file x2, .stations file

% Aim of this script: Go from a coordinate-based file to an IGMAS modelling
% file

% Advantages: intuitive to model, all possible caves can be modelled,
% script is relatively understandable (apart from some technicalities)
% Disadvantage: making the grid takes quite some time, voxels must be
% placed on a grid (so keep enough distance between voxels), when loaded,
% one can not check if the cave is properly modelled

% What can be changed?
% - Reducing factor: how many data points do you find acceptable (Step 0)
% - Cavity depth: The distance between cavity and stations (Step 1)
% - Rock density: what the density of the surrounding rock is (Step 1)
% - Sizing factor: How much larger do you want the model to be, compared to
% the cavity? To deal with side effects, but also increases complexity
% (Step 1)
% - Station resolution: how many stations do you want on the grid (Step 1)
% - Voxel resolution: how many voxels do you maximally (!) want per
% direction? It is made such that the grid has equal spacing in all
% directions (Step 1)
% - Name of the model (Step 1)

clear all
close all
%% Step 0: Import the coordinates and reduce the resolution

% the file either is an csv-file (use importdata) or a .stl file (use stlread)

% Cavity_coordinates = importdata('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\MatLab Model Maker\Cavity_Test_Data_Coords_4.csv'); % Change the directory accordingly
[A,~] = stlread('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Lava tube models\galapagos-lava-tube\source\Mesh SCruz2 levigata\Mesh SCruz2 levigata.stl');
Cavity_coordinates = A.Points;
Reduce_factor = 5;

% What this does; make the cavity Sizing_Factor times larger, or stretch it
% by a factor Stretching_Factor_H horizontally.
% Vertical stretching is a bit trickier, or shape-preserving stretching
% Alternative; stretch it until the average radius/length is .../... meter
Sizing_Factor = 1;
Stretching_Factor_H = 1;

Cavity_coordinates = Sizing_Factor*Reduce_coordinates(Cavity_coordinates, Reduce_factor);
Cavity_coordinates(:,1) = Stretching_Factor_H*Cavity_coordinates(:,1);
Cavity_coordinates(:,2) = Stretching_Factor_H*Cavity_coordinates(:,2);

%% Step 1: Export a basic .model file
% Main properties: dimensions of whole model, bodies

Rock_density = 2.5; % density in g/cm³
Cavity_depth = 10; % How many meters under the surface the cave is located

Station_resolution = 50; % Square root how number of stations wanted
Voxel_resolution = 100; % Number of voxels per direction. Note: in total 30.000.000 voxels are allowed.

Sizing_factor = 1.50; % To deal with edge effects

Model_name = ['New_voxle_test']; % The name you want your model to have

Cavity_limit_vector = [];
for a=1:3
    Cavity_limit_vector = [Cavity_limit_vector,round(min(Cavity_coordinates(:,a))),round(max(Cavity_coordinates(:,a)))];
end

% Cavity_length_x = Cavity_limit_vector(2)-Cavity_limit_vector(1);
% Cavity_length_y = Cavity_limit_vector(4)-Cavity_limit_vector(3);

Cavity_coordinates(:,3) = Cavity_coordinates(:,3) - Cavity_limit_vector(6) - Cavity_depth;

Cavity_limit_vector(5) = Cavity_limit_vector(5) - Cavity_limit_vector(6)-Cavity_depth;
Cavity_limit_vector(6) = 0;

Meter_per_voxel_x = (Cavity_limit_vector(2)-Cavity_limit_vector(1))/Voxel_resolution;
Meter_per_voxel_y = (Cavity_limit_vector(4)-Cavity_limit_vector(3))/Voxel_resolution;
Meter_per_voxel_z = (Cavity_limit_vector(6)-Cavity_limit_vector(5))/Voxel_resolution;
Grid_spacing = round(max([Meter_per_voxel_x,Meter_per_voxel_y,Meter_per_voxel_z]));

Cavity_limit_vector(1) = Cavity_limit_vector(1) - mod((Cavity_limit_vector(2)-Cavity_limit_vector(1)),Grid_spacing);
Cavity_limit_vector(3) = Cavity_limit_vector(3) - mod((Cavity_limit_vector(4)-Cavity_limit_vector(3)),Grid_spacing);
Cavity_limit_vector(5) = Cavity_limit_vector(5) - mod((Cavity_limit_vector(6)-Cavity_limit_vector(5)),Grid_spacing);

% Now, everything fits into the grid

Voxels_in_direction_x = (Cavity_limit_vector(2)-Cavity_limit_vector(1))/Grid_spacing;
Voxels_in_direction_y = (Cavity_limit_vector(4)-Cavity_limit_vector(3))/Grid_spacing;
Voxels_in_direction_z = (Cavity_limit_vector(6)-Cavity_limit_vector(5))/Grid_spacing;

Number_voxels_per_direction = [Voxels_in_direction_x,Voxels_in_direction_y,Voxels_in_direction_z];

%% Step 2: Create a voxel file for the cave, and one without any cave for comparison
Vox_limits = [];
% With cavity
File_name = ['Vox_With_Cavity_',Model_name,'.vxo'];
Cavity_bool = true;
Vox_limits = Voxel_writing(Rock_density,Cavity_coordinates,Cavity_limit_vector,Cavity_bool,File_name,Grid_spacing,Number_voxels_per_direction,Vox_limits,Sizing_factor,Model_name);

% Without cavity
File_name = ['Vox_No_Cavity_',Model_name,'.vxo'];
Cavity_bool = false;
Vox_limits = Voxel_writing(Rock_density,Cavity_coordinates,Cavity_limit_vector,Cavity_bool,File_name,Grid_spacing,Number_voxels_per_direction,Vox_limits,Sizing_factor,Model_name);

%% Step 3: Define the stations

Export_stations_file(Cavity_limit_vector,Station_resolution,Model_name) % ASSUMING flat land, also currently only a set grid is created without other patterns

%% Step 4: Import the model file, the voxel files and the station file in IGMAS, do the anomaly calculation and export the .stations files

%% Step 5: Compare the calculated anomalies

% Import the .stations file with results. The 'xml2struct' function is a
% function from file exchange
% (https://nl.mathworks.com/matlabcentral/fileexchange/28518-xml2struct)
% which makes it easy to extract data (coordinates, values) from files like
% the .stations file.

% Put this at true if you have results to analyse
Calculation_bool = false;

if Calculation_bool == true
    With_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package\Sizing Tests and Results\f10_With_Cavity.stations');

    With_cavity_results_array = [];
    for a=1:length(With_cavity_results.geodata.vertex(1,:))

        Temp_x = With_cavity_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = With_cavity_results.geodata.vertex{1,a}.Attributes.y;

        Temp_result = With_cavity_results.geodata.vertex{1,a}.property.Attributes.value;

        With_cavity_results_array(a,1) = str2num(Temp_x);
        With_cavity_results_array(a,2) = str2num(Temp_y);
        With_cavity_results_array(a,3) = str2num(Temp_result);
    end
    clear a

    No_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package\Sizing Tests and Results\f10_No_Cavity.stations');

    No_cavity_results_array = [];
    for a=1:length(No_cavity_results.geodata.vertex(1,:))

        Temp_x = No_cavity_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = No_cavity_results.geodata.vertex{1,a}.Attributes.y;

        Temp_result = No_cavity_results.geodata.vertex{1,a}.property.Attributes.value;

        No_cavity_results_array(a,1) = str2num(Temp_x);
        No_cavity_results_array(a,2) = str2num(Temp_y);
        No_cavity_results_array(a,3) = str2num(Temp_result);
    end
    clear a

    % With the results of both files, we can make a graph of the difference
    % between the models with and without cavity. Note that the range is larger
    % than the one of the cavity for reducing edge effects. However, since the
    % edge effect is present in both models, taking the difference (almost) eliminates
    % that.
    
%     figure()
%     scatter3(With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3))
%     title('Difference between the signals of models with and without a cavity','interpreter','latex')
%     xlabel('x-coordinate [m]','interpreter','latex')
%     ylabel('y-coordinate [m]','interpreter','latex')
%     zlabel('Gravity Difference [mGal]','interpreter','latex')
%     
%     figure()
%     stem3(With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3))
%     title('Difference between the signals of models with and without a cavity','interpreter','latex')
%     xlabel('x-coordinate [m]','interpreter','latex')
%     ylabel('y-coordinate [m]','interpreter','latex')
%     zlabel('Gravity Difference [mGal]','interpreter','latex')
    
    tri = delaunay(With_cavity_results_array(:,1), With_cavity_results_array(:,2));
    figure()
    axis tight
    trisurf(tri, With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3))
    title('Difference between the signals of models with and without a cavity, Sizing 10','interpreter','latex')
    xlabel('x-coordinate [m]','interpreter','latex')
    ylabel('y-coordinate [m]','interpreter','latex')
    zlabel('Gravity Difference [mGal]','interpreter','latex')
    c = colorbar;
    c.Label.String = 'Gravity Difference [mGal]';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 11;
    
    labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [mGal]','Difference between the signals of models with and without a cavity, Sizing 10' };
    ContourPlot(2,With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3),labels)
    
    
    
    With_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package\Sizing Tests and Results\f1_With_Cavity.stations');

    With_cavity_results_array = [];
    for a=1:length(With_cavity_results.geodata.vertex(1,:))

        Temp_x = With_cavity_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = With_cavity_results.geodata.vertex{1,a}.Attributes.y;

        Temp_result = With_cavity_results.geodata.vertex{1,a}.property.Attributes.value;

        With_cavity_results_array(a,1) = str2num(Temp_x);
        With_cavity_results_array(a,2) = str2num(Temp_y);
        With_cavity_results_array(a,3) = str2num(Temp_result);
    end
    clear a

    No_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package\Sizing Tests and Results\f1_No_Cavity.stations');

    No_cavity_results_array = [];
    for a=1:length(No_cavity_results.geodata.vertex(1,:))

        Temp_x = No_cavity_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = No_cavity_results.geodata.vertex{1,a}.Attributes.y;

        Temp_result = No_cavity_results.geodata.vertex{1,a}.property.Attributes.value;

        No_cavity_results_array(a,1) = str2num(Temp_x);
        No_cavity_results_array(a,2) = str2num(Temp_y);
        No_cavity_results_array(a,3) = str2num(Temp_result);
    end
    clear a

    % With the results of both files, we can make a graph of the difference
    % between the models with and without cavity. Note that the range is larger
    % than the one of the cavity for reducing edge effects. However, since the
    % edge effect is present in both models, taking the difference (almost) eliminates
    % that.
    
%     figure()
%     scatter3(With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3))
%     title('Difference between the signals of models with and without a cavity','interpreter','latex')
%     xlabel('x-coordinate [m]','interpreter','latex')
%     ylabel('y-coordinate [m]','interpreter','latex')
%     zlabel('Gravity Difference [mGal]','interpreter','latex')
%     
%     figure()
%     stem3(With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3))
%     title('Difference between the signals of models with and without a cavity','interpreter','latex')
%     xlabel('x-coordinate [m]','interpreter','latex')
%     ylabel('y-coordinate [m]','interpreter','latex')
%     zlabel('Gravity Difference [mGal]','interpreter','latex')
    
    tri = delaunay(With_cavity_results_array(:,1), With_cavity_results_array(:,2));
    figure()
    axis tight
    trisurf(tri, With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3))
    title('Difference between the signals of models with and without a cavity, Sizing 1','interpreter','latex')
    xlabel('x-coordinate [m]','interpreter','latex')
    ylabel('y-coordinate [m]','interpreter','latex')
    zlabel('Gravity Difference [mGal]','interpreter','latex')
    c = colorbar;
    c.Label.String = 'Gravity Difference [mGal]';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 11;
    
    labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [mGal]','Difference between the signals of models with and without a cavity, Sizing 1' };
    ContourPlot(4,With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3),labels)
end

%% Functions

% Function to also export a stations file. There is no elevation assumed,
% and also just a standard grid. This can be improved in the future.
function Export_stations_file(Station_limit_vector,Station_resolution, Model_name)
    feature('DefaultCharacterSet', 'UTF8'); % Needed mostly to make the ³-character work...
    File_name = ['Vox_Stations_',Model_name,'.stations'];
    fid = fopen(File_name,'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n'); % I am not sure what this does; result in proper encoding?
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

%% functions
function ContourPlot(fignumber,x,y,z,labels)
figure(fignumber)
xlin = linspace(min(x),max(x),2500);                         
ylin = linspace(min(y),max(y),2500);                         
[X,Y] = meshgrid(xlin,ylin);                                      
Z = griddata(x,y,z,X,Y,'linear');
contourf(X,Y,Z, 'LineColor', 'none' )
hold on
xlabel(labels(1),'interpreter','latex','fontsize',12)
ylabel(labels(2),'interpreter','latex','fontsize',12)
cb=colorbar;
cb.Label.String = labels(3);
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 12;
caxis([min(z) max(z)])
title(labels(4),'interpreter','latex','fontsize',12)
set(gcf,'position',[100,100,800,400])
end

function ContourSubPlot(fignumber,x,y,z,labels)
figure(fignumber)
xlin = linspace(min(x),max(x),2500);                         
ylin = linspace(min(y),max(y),2500);                         
[X,Y] = meshgrid(xlin,ylin);                                      
Z = griddata(x,y,z,X,Y,'linear');
contourf(X,Y,Z, 'LineColor', 'none' )
hold on
xlabel(labels(1),'interpreter','latex','fontsize',12)
ylabel(labels(2),'interpreter','latex','fontsize',12)
cb=colorbar;
cb.Label.String = labels(3);
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 12;
caxis([min(z) max(z)])
set(gcf,'position',[100,100,600,300])
end


