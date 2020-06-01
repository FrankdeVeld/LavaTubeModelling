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

% Cavity_coordinates = importdata('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\MatLab Model Maker\Cavity_Test_Data_Coords_4.csv'); % Change the directory accordingly

[A,~] = stlread('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Lava tube models\cueva-de-los-siete-lagos\source\agua_text\Los lagos ok.stl');
Cavity_coordinates = A.Points;
Reduced_index = 1;
Reduce_factor = 50;
Reduced_cavity_coordinates = zeros(floor(length(Cavity_coordinates(:,1))/Reduce_factor),3);
for a = 1:length(Cavity_coordinates(:,1))
    if mod(a,Reduce_factor)==0
        Reduced_cavity_coordinates(Reduced_index,:) = Cavity_coordinates(a,:);
        Reduced_index = Reduced_index + 1;
    end
end
clear a
Cavity_coordinates = Reduced_cavity_coordinates;

%% Step 1: Export a basic .model file
% Main properties: dimensions of whole model, bodies

Rock_density = 2.5; % density in g/cm³
Cavity_depth = 10; % How many meters under the surface the cave is located
Sizing_factor = 1.5; % To deal with edge effects
Station_resolution = 50; % Square root how number of stations wanted
Voxel_resolution = 500; % Number of voxels per direction. Note: in total 30.000.000 voxels are allowed.

Model_name = ['Comparison_Test']; % The name you want your model to have

Cavity_limit_vector = [];
for a=1:3
    Cavity_limit_vector = [Cavity_limit_vector,min(Cavity_coordinates(:,a)),max(Cavity_coordinates(:,a))];
end
% Interpretation: the total model has the size of the cavity, times 2*Sinzing_factor.
% This can be tweaked!
Cavity_length_x = Cavity_limit_vector(2)-Cavity_limit_vector(1);
Cavity_length_y = Cavity_limit_vector(4)-Cavity_limit_vector(3);

Limit_vector(1) = round(Cavity_limit_vector(1) - Sizing_factor*Cavity_length_x);
Limit_vector(2) = round(Cavity_limit_vector(2) + Sizing_factor*Cavity_length_x);
Limit_vector(3) = round(Cavity_limit_vector(3) - Sizing_factor*Cavity_length_y);
Limit_vector(4) = round(Cavity_limit_vector(4) + Sizing_factor*Cavity_length_y);
Limit_vector(5) = round(Cavity_limit_vector(5) - Cavity_depth); % We need a 'buffer' here
Limit_vector(6) = round(Cavity_limit_vector(6) + Cavity_depth);

Cavity_coordinates(:,3) = Cavity_coordinates(:,3) - Limit_vector(6);
Limit_vector(5)  = Limit_vector(5) - Limit_vector(6);
Limit_vector(6)  = Limit_vector(6) - Limit_vector(6);

res_x = (Limit_vector(2)-Limit_vector(1))/Voxel_resolution;
res_y = (Limit_vector(4)-Limit_vector(3))/Voxel_resolution;
res_z = (Limit_vector(6)-Limit_vector(5))/Voxel_resolution;
Grid_size = round(max([res_x,res_y,res_z]));

Limit_vector(1) = Limit_vector(1) + mod((Limit_vector(2)-Limit_vector(1)),Grid_size);
Limit_vector(3) = Limit_vector(3) + mod((Limit_vector(4)-Limit_vector(3)),Grid_size);
Limit_vector(5) = Limit_vector(5) + mod((Limit_vector(6)-Limit_vector(5)),Grid_size);

% Now, everything fits into the grid

res_x = (Limit_vector(2)-Limit_vector(1))/Grid_size;
res_y = (Limit_vector(4)-Limit_vector(3))/Grid_size;
res_z = (Limit_vector(6)-Limit_vector(5))/Grid_size;

Grid = [res_x,res_y,res_z];

Write_basic_model_file(Rock_density, Model_name, Limit_vector, res_y-2);
% Note: this model has almost nothing; the proper cavity will be defined
% later with voxels

%% Step 2: Create a voxel file for the cave, and one without any cave for comparison

% With cavity
File_name = ['Vox_Voxel_with_cavity_',Model_name,'.vxo'];
Cavity_bool = true;
Voxel_writing(Voxel_resolution,Rock_density,Cavity_coordinates,Limit_vector,Cavity_bool,File_name,Grid_size,Grid)

% Without cavity
File_name = ['Vox_Voxel_no_cavity_',Model_name,'.vxo'];
Cavity_bool = false;
Voxel_writing(Voxel_resolution,Rock_density,Cavity_coordinates,Limit_vector,Cavity_bool,File_name,Grid_size,Grid)

%% Step 3: Define the stations
% In the same way as in the definition of the model dimensions, here the
% station grid location is defined. However, it turned out to be better to
% measure everywhere
% Station_limit_vector(1) = Cavity_limit_vector(1) - Station_factor*Cavity_length_x;
% Station_limit_vector(2) = Cavity_limit_vector(2) + Station_factor*Cavity_length_x;
% Station_limit_vector(3) = Cavity_limit_vector(3) - Station_factor*Cavity_length_y;
% Station_limit_vector(4) = Cavity_limit_vector(4) + Station_factor*Cavity_length_y;
% Station_limit_vector(5) = Limit_vector(5);
% Station_limit_vector(6) = Limit_vector(6);
% Idea: you need to measure even at the boundaries, and later you limit
% this to the above limits

Export_stations_file(Limit_vector,Station_resolution,Model_name) % ASSUMING flat land, also currently only a set grid is created without other patterns

%% Step 4: Import the model file, the voxel files ánd the station file in IGMAS, do the anomaly calculation and export the .stations files

%% Step 5: Compare the calculated anomalies

% Import the .stations file with results. The 'xml2struct' function is a
% function from file exchange
% (https://nl.mathworks.com/matlabcentral/fileexchange/28518-xml2struct)
% which makes it easy to extract data (coordinates, values) from files like
% the .stations file.

% Put this at true if you ahve results to analyse
Calculation_bool = false;

if Calculation_bool == true
    With_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\TestsMay27\Vox_With_cavity_results_Los_Lagos_May28_reduce50_vox250.stations');

    With_cavity_results_array = [];
    for a=1:length(With_cavity_results.geodata.vertex(1,:))

        Temp_x = With_cavity_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = With_cavity_results.geodata.vertex{1,a}.Attributes.y;

        Temp_result = With_cavity_results.geodata.vertex{1,a}.property{1,1}.Attributes.value;

        With_cavity_results_array(a,1) = str2num(Temp_x);
        With_cavity_results_array(a,2) = str2num(Temp_y);
        With_cavity_results_array(a,3) = str2num(Temp_result);
    end
    clear a

    No_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\TestsMay27\Vox_No_cavity_results_Los_Lagos_May28_reduce50_vox250.stations');

    No_cavity_results_array = [];
    for a=1:length(No_cavity_results.geodata.vertex(1,:))

        Temp_x = No_cavity_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = No_cavity_results.geodata.vertex{1,a}.Attributes.y;

        Temp_result = No_cavity_results.geodata.vertex{1,a}.property{1,1}.Attributes.value;

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
    figure()
    scatter3(With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3))
    title('Difference between the signals of models with and without a cavity','interpreter','latex')
    xlabel('x-coordinate [m]','interpreter','latex')
    ylabel('y-coordinate [m]','interpreter','latex')
    zlabel('Gravity Difference [mGal]','interpreter','latex')

end

%% Aiding functions
% Function to write a basic .model file, according to the previously
% defined model limits
function Write_basic_model_file(Rock_density, Model_name, Limit_vector, Num_intermediate_sections) 
    feature('DefaultCharacterSet', 'UTF8'); % Needed mostly to make the ³-character work
    File_name = ['Vox_Basic_model_',Model_name,'.model'];
    fid = fopen(File_name,'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n'); % I am not sure what this does; result in proper encoding?
    fprintf(fid, '<geodata name="Test Model">\n'); % Open the geodata class
    fprintf(fid, '<projection name="unknown" units="m"></projection> \n');
 
    
    % Define a reference body
    fprintf(fid, '<property name="body" value="reference"> \n    <property name="density" units="g/cm³" value="0.0"></property> \n <color red="0.5019608" green="0.5019608" blue="0.5019608"></color> \n</property> \n' );

    % Define another, rocky body. Note: density can be changed!
    String_rock = ['<property name="body" value="stone"> \n    <property name="density" units="g/cm³" value="%4.2f"></property> \n <color red="1" green="0" blue="0"></color> \n</property> \n'];
    String_rock = sprintf(String_rock,Rock_density);
    fprintf(fid, String_rock);

    % Define the cavity body. Density should be 0.
    String_rock = ['<property name="body" value="cavity"> \n    <property name="density" units="g/cm³" value="0"></property> \n <color red="0" green="1" blue="0"></color> \n</property> \n'];
    fprintf(fid, String_rock);
    
    % Section 0
    
    Section_text = ['<geometry>' ...
      '<cross_section name="0" '...
      'x_start="%4.2f" y_start="%4.2f" x_end="%4.2f" y_end="%4.2f"></cross_section> \n'];
    Section_text = sprintf(Section_text,Limit_vector(1),Limit_vector(3),Limit_vector(2),Limit_vector(3));

    Temp_text = ['<vertex id="0" x="%4.2f" z="%4.2f"></vertex> \n' ...
    '<vertex id="1" x="%4.2f" z="%4.2f"></vertex> \n' ....
    '<vertex id="2" x="%4.2f" z="%4.2f"></vertex> \n'...
    '<vertex id="3" x="%4.2f" z="%4.2f"></vertex> \n'];
    Temp_text = sprintf(Temp_text,Limit_vector(2)-Limit_vector(1),Limit_vector(5),0,Limit_vector(5),0,Limit_vector(6),Limit_vector(2)-Limit_vector(1),Limit_vector(6));
    Section_text = [Section_text,Temp_text];

    Section_text = [Section_text,'<entry type="polygon" id_list="0 1 2 3 "> \n' ...
         '<property name="body" value="stone"></property> \n'...
      '</entry> \n']; 

    Section_text = [Section_text,'</geometry> \n'];
    fprintf(fid,Section_text);
    
    index = 3;
    
    for a=1:Num_intermediate_sections
        Section_text = ['<geometry>' ...
          '<cross_section name="%d" '...
          'x_start="%4.2f" y_start="%4.2f" x_end="%4.2f" y_end="%4.2f"></cross_section> \n'];
        Section_text = sprintf(Section_text,a,Limit_vector(1),a*(Limit_vector(4)-Limit_vector(3))/(Num_intermediate_sections + 1)+ Limit_vector(3),Limit_vector(2),a*(Limit_vector(4)-Limit_vector(3))/(Num_intermediate_sections + 1)+ Limit_vector(3));

        Temp_text = ['<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n' ...
        '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n' ....
        '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'...
        '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'];
        Temp_text = sprintf(Temp_text,index + 1, Limit_vector(2)-Limit_vector(1),Limit_vector(5),index + 2,0,Limit_vector(5),index + 3,0,Limit_vector(6),index + 4,Limit_vector(2)-Limit_vector(1),Limit_vector(6));
        Section_text = [Section_text,Temp_text];
        
        Temp_text = ['<entry type="polygon" id_list="%d %d %d %d "> \n' ...
             '<property name="body" value="stone"></property> \n'...
          '</entry> \n'];
        Temp_text = sprintf(Temp_text,index+1,index+2,index+3,index+4);
        Section_text = [Section_text,Temp_text];
        
        index = index + 4;
        
        Section_text = [Section_text,'</geometry> \n'];
        fprintf(fid,Section_text);
    end
    
    % Last section
    
    Section_text = ['<geometry>' ...
      '<cross_section name="%d" '...
      'x_start="%4.2f" y_start="%4.2f" x_end="%4.2f" y_end="%4.2f"></cross_section> \n'];
    Section_text = sprintf(Section_text,a+1,Limit_vector(1),Limit_vector(4),Limit_vector(2),Limit_vector(4));

    Temp_text = ['<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n' ...
    '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n' ....
    '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'...
    '<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'];
    Temp_text = sprintf(Temp_text,index + 1, Limit_vector(2)-Limit_vector(1),Limit_vector(5),index + 2, 0,Limit_vector(5),index + 3,0,Limit_vector(6),index + 4,Limit_vector(2)-Limit_vector(1),Limit_vector(6));
    Section_text = [Section_text,Temp_text];
    
    Temp_text = ['<entry type="polygon" id_list="%d %d %d %d "> \n' ...
         '<property name="body" value="stone"></property> \n'...
      '</entry> \n'];
    Temp_text = sprintf(Temp_text,index+1,index+2,index+3,index+4);
    Section_text = [Section_text, Temp_text];

    Section_text = [Section_text,'</geometry> \n'];
    fprintf(fid,Section_text);


    fprintf(fid, '</geodata>\n' );  % Close the geodata class, and with this the file
    fclose(fid);
end

% Function to write .vxo-files
function  Voxel_writing(Voxel_resolution,Rock_density,Cavity_coordinates,Limit_vector,Cavity_bool,File_name,Grid_size,Grid)
    Voxel = []; % Initialize the full data file
    xlength = Voxel_resolution; % Voxel_resolution: number of voxels per direction
    ylength = Voxel_resolution; 
    zlength = Voxel_resolution; 
    
    xrange = Limit_vector(2) - Limit_vector(1);
    yrange = Limit_vector(4) - Limit_vector(3);
    zrange = Limit_vector(6) - Limit_vector(5);
    index = 1;
    Count_index = 0;
    
    Direct_grid_bool = false;
    Middle_search_bool = false;
    Smart_grid_bool = true;
    
    % Method 1: gridded network of voxels. Takes too long, too much rock.
    
    if Direct_grid_bool==true
        for x=0:xlength 
            for y=0:ylength 
                for z=0:zlength 
                    Voxel(index,1) = Limit_vector(1) + xrange/Voxel_resolution*x + 0.5; % x-coordinate 
                    Voxel(index,2) = Limit_vector(3) + yrange/Voxel_resolution*y + 0.5; % y-coordinate
                    Voxel(index,3) = Limit_vector(5) + zrange/Voxel_resolution*z + 0.5; % z-coordinate

                    Current_coordinate = [Voxel(index,1),Voxel(index,2),Voxel(index,3)];
                    [Voxel(index,4),Count_index] = Calculate_density(Rock_density,Cavity_coordinates,Cavity_bool,Current_coordinate,Count_index); % Density
                    
                    index = index + 1;
                end
            end
        end
        
        [Voxel,index] = Add_corner_voxels(index,Voxel,Limit_vector,Rock_density);
    end
    
    % Method 2: make 'fake sections', round coordinates to these sections,
    % find the middle points and start from there
    
    % Issue: IGMAS wants voxels on grids, apparently. Thus, we round all of
    % these to a grid
    if Middle_search_bool == true
        Number_cavity_sections = round((max(Cavity_coordinates(:,2))-min(Cavity_coordinates(:,2)))/Grid_size);
        Number_cavity_sections = Number_cavity_sections - mod(Number_cavity_sections,Grid_size);
        Number_samples_per_section = 50;
        Sorted_cavity_coordinates = Fake_sectioning(Number_cavity_sections,Cavity_coordinates,yrange);
        
        Coordinate_cell = Reshaper(Number_cavity_sections,Sorted_cavity_coordinates);
        
        
        
        % Determines the ratio rock voxels to cavity voxels
        Sampling_radius_factor = 1;
        
        Middle_points = [];
        for a=1:Number_cavity_sections
            [x_middle, z_middle] = Middle_point_seeker(cell2mat(Coordinate_cell(a,1)));
            y_middle = (a-1)*Number_cavity_sections/yrange;
            New_middle_point = [x_middle, y_middle, z_middle];
            Middle_points = [Middle_points, New_middle_point];
            
            Max_distance = Distance_calculator(New_middle_point,cell2mat(Coordinate_cell(a,1)));
            Sampling_radius = Sampling_radius_factor * Max_distance;
            for b=1:Number_samples_per_section
                % Generate random points on a circle
                theta_sample = rand*2*pi;
                r_sample = rand*Sampling_radius^2;
                
                x_sample = sqrt(r_sample)*cos(theta_sample);
                z_sample = sqrt(r_sample)*sin(theta_sample);
                
                y_sample = y_middle;
                
%                 if mod(xrange,x_sample) < 0.5*Grid_size
%                     x_sample = x_sample - mod(xrange,x_sample);
%                 else
%                     x_sample = x_sample + mod(xrange,x_sample);
%                 end
%                 
%                 if mod(yrange,y_sample) < 0.5*Grid_size
%                     y_sample = y_sample - mod(xrange,y_sample);
%                 else
%                     y_sample = y_sample + mod(xrange,y_sample);
%                 end
%                 
%                 if mod(yrange,z_sample) < 0.5*Grid_size
%                     z_sample = z_sample - mod(xrange,z_sample);
%                 else
%                     z_sample = z_sample + mod(xrange,z_sample);
%                 end
                
                % Now the newly generated coordinates are also on the grid
                
                Voxel(index,1) = x_sample;
                Voxel(index,2) = y_sample;
                Voxel(index,3) = z_sample;
                
                Current_coordinate = [Voxel(index,1),Voxel(index,2),Voxel(index,3)];
                [Voxel(index,4),Count_index] = Calculate_density(Rock_density,Cavity_coordinates,Cavity_bool,Current_coordinate,Count_index); % Density
                
                index = index + 1;
            end
            clear b
        end
        clear a
        
        [Voxel,index] = Add_corner_voxels(index,Voxel,Limit_vector,Rock_density);
    end
     
    % Method 3: make a grid of rock, locally probe the cavity
    if Smart_grid_bool == true
        Number_samples_per_section = 5;
        
        Starting_grid = Make_rocky_grid(Limit_vector,Rock_density,Grid_size,Grid);

        Reduced_cavity_coordinates = Grid_cavity_coordinates(Cavity_coordinates,Grid_size,Starting_grid);

        Unique_reduced_cavity_coordinates = unique(Reduced_cavity_coordinates, 'rows');

        [~,idx] = sort(Unique_reduced_cavity_coordinates(:,2)); % sort just the first column
        Sorted_unique_reduced_cavity_coordinates = Unique_reduced_cavity_coordinates(idx,:);   % sort the whole matrix using the sort indices
        
        % Due to technical issues, we rather don't have zeros in the
        % matrix. If this is really a problem, a method can be implemented
        % that temporarily changes the zero coordinates to another value.
        for a=1:length(Sorted_unique_reduced_cavity_coordinates(:,1))
            if Sorted_unique_reduced_cavity_coordinates(a,1) == 0
                Sorted_unique_reduced_cavity_coordinates(a,1) = Sorted_unique_reduced_cavity_coordinates(a,1) + 0.001;
            end
            if Sorted_unique_reduced_cavity_coordinates(a,2) == 0
                Sorted_unique_reduced_cavity_coordinates(a,2) = Sorted_unique_reduced_cavity_coordinates(a,2) + 0.001;
            end
            if Sorted_unique_reduced_cavity_coordinates(a,3) == 0
                Sorted_unique_reduced_cavity_coordinates(a,3) = Sorted_unique_reduced_cavity_coordinates(a,3) + 0.001;
            end
        end
        
        Sorted_reduced_cavity_y_coordinates = Sorted_unique_reduced_cavity_coordinates(:,2);
        Sorted_unique_reduced_cavity_y_coordinates = unique(Sorted_reduced_cavity_y_coordinates);
        Number_cavity_sections = length(Sorted_unique_reduced_cavity_y_coordinates);
        Coordinate_cell = Reshaper(Number_cavity_sections,Sorted_unique_reduced_cavity_coordinates);
        

        
        Middle_points = [];
        for a=1:Number_cavity_sections
            Current_section_matrix = cell2mat(Coordinate_cell(a,1));
            
            % As it turns out, having values like 0.001 in the matrix indeed is
            % a problem. We transfer them back to zeros
            for b=1:length(Current_section_matrix(:,1))
                if Current_section_matrix(b,1) == 0.001
                    Current_section_matrix(b,1) = 0;
                end
                if Current_section_matrix(b,2) == 0.001
                    Current_section_matrix(b,2) = 0;
                end
            end
            clear b
            
            Max_x_distance = Grid_distance_calculator(Current_section_matrix(:,1));
            Max_z_distance = Grid_distance_calculator(Current_section_matrix(:,2));
            
            Number_x_point_evaluation = Max_x_distance/Grid_size + 1;
            Number_z_point_evaluation = Max_z_distance/Grid_size + 1;
            
            Min_x = min(Current_section_matrix(:,1));
            Min_z = min(Current_section_matrix(:,2));
            
            for b=1:Number_x_point_evaluation
                for c=1:Number_z_point_evaluation
                    Temp_x = Min_x + b*Grid_size;
                    Temp_z = Min_z + c*Grid_size;
                    Temp_y = Sorted_unique_reduced_cavity_y_coordinates(a);
                    
                    if Temp_y == 0.001
                        Temp_y = 0;
                    end
                    
                    Temp_coord = [Temp_x,Temp_y,Temp_z]; 
                    
                    [Temp_density,Count_index] = Calculate_density(Rock_density,Cavity_coordinates,Cavity_bool,Temp_coord,Count_index); % Density
                    
                    Temp_coord = [Temp_coord,2.5];
                    if Temp_density == 0
                        [~,~,index] = intersect(Temp_coord,Starting_grid,'rows');
                        Temp_coord = [Temp_x,Temp_y,Temp_z,Temp_density];
                        for e =1:4
                            Starting_grid(index,e) = Temp_coord(e);
                        end
                    end
                end 
            end
        end
        Voxel = Starting_grid;
    end
    disp(Count_index)
    

    % Writing it to a .vxo file
    fileID = fopen(File_name,'w');
    fprintf(fileID,'x\t y\t z\t cellValue\n');
    for i = 1:length(Voxel(:,1))
      fprintf(fileID,'%d\t%d\t%d\t%d\n', Voxel(i,1),Voxel(i,2),Voxel(i,3),Voxel(i,4));
    end
    fclose(fileID);
end

% Function to write a density to the voxels; either Rock_density or 0
function [Density,Count_index] = Calculate_density(Rock_density,Cavity_coordinates,Cavity_bool,Current_coordinate,Count_index)
    if Cavity_bool == false
        Density = Rock_density;
    else
        % One of many ways to check whether a point is inside a polyhedron.
        % Sample code obtained from https://nl.mathworks.com/matlabcentral/answers/101396-is-there-a-function-in-matlab-for-detecting-points-inside-a-polyhedron
        % 'Delaunayn' and 'Tsearchn' are MatLab standard funcions
        Triangularization = delaunayn([Cavity_coordinates(:,1) Cavity_coordinates(:,2) Cavity_coordinates(:,3)]); % Generate delaunay triangulization
        Search = tsearchn([Cavity_coordinates(:,1) Cavity_coordinates(:,2) Cavity_coordinates(:,3)], Triangularization, Current_coordinate); % Determine which triangle point is within
        Inside_bool = ~isnan(Search); % Convert to logical vector
        
        if Inside_bool == true
            Density = 0;
            Count_index = Count_index + 1;
        else
            Density = Rock_density;
        end
    end
    % Maybe interesting alternative:
    % https://nl.mathworks.com/matlabcentral/answers/152189-volume-of-3d-polyhedron
end

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

% Function to re-arrange the coordinates like the sectioning program
function Sorted_cavity_coordinates = Fake_sectioning(Number_fake_sections,Cavity_coordinates,y_range)
    Leftover = mod(length(Cavity_coordinates(:,2)),Number_fake_sections); % Idea: assign x points to one section, but with this you have leftovers
    Leftover_front = ceil(Leftover/2); % Crude method: put half of the left-overs in the first section
    Leftover_back = floor(Leftover/2); % Crude method: put half of the left-overs in the last section

    [~,Sorted_y] = sort(Cavity_coordinates(:,2)); % Sort the y-coordinates
    Sorted_cavity_coordinates = Cavity_coordinates(Sorted_y,:); % Sort the whole matrix according to these y-coordinates

    % It is chosen to have ymin, ymax and (Number_of_sections - 1)
    % points defining bins, so that the middle of each bin is defined
    % by a section
    % Like this:
    % Number_of_sections = 10
    % |.....|.....|.....|.....|.....|.....|.....|.....|.....|.....|
    % ymin .1.....2.....3.....4.....5.....6.....7.....8.....9..ymax
    % Sections:
    % ...^.....^.....^.....^.....^.....^.....^.....^.....^.....^...
    % ..s1....s2....s3....s4....s5....s6....s7....s8....s9....s10..
    for a=1:Leftover_front
        y_coor = min(Cavity_coordinates(:,2)) + 0.5*y_range/Number_fake_sections; % Thus, the leftovers at the front are combined in the first section
        Sorted_cavity_coordinates(a,2) = y_coor;
    end
    clear a

    for a=(Leftover_front+1):(length(Sorted_cavity_coordinates(:,2))-Leftover_back-1) % This is a loop through the sorted coordinates
        Nodes_per_section = length(Cavity_coordinates(:,1))/Number_fake_sections;
        Current_section = floor(a/Nodes_per_section); % The information on the section we are working in is still needed
        y_coor = min(Cavity_coordinates(:,2)) + y_range/Number_fake_sections*Current_section + 0.5*y_range/Number_fake_sections;
        % Explanation: this is the y-value, so the minimum is
        % min(Cavity_coordinates(:,2)) and the middle of the current
        % section is half a section away from the 'nodes' explained
        % before. Lastly, y_range/Number_of_sections*Current_section is
        % to account for previous sections
        Sorted_cavity_coordinates(a,2) = y_coor;
    end
    clear a

    for a=(length(Sorted_cavity_coordinates(:,2))-Leftover_back):length(Sorted_cavity_coordinates(:,2)) % The same thing happens for the back as happened for the front
        y_coor = min(Cavity_coordinates(:,2)) + y_range - 0.5*y_range/Number_fake_sections;
        Sorted_cavity_coordinates(a,2) = y_coor;
    end
    clear a

end

% Function to find the middle of a cavity in a 'fake' section
function [Av_x, Av_z] = Middle_point_seeker(Section_coords) 
    Tot_x = 0;
    Tot_z = 0;
    for a=1:length(Section_coords(:,1))
        Tot_x = Tot_x + Section_coords(a,1);
        Tot_z = Tot_z + Section_coords(a,2);
    end
    clear a
    Av_x = Tot_x/length(Section_coords(:,1));
    Av_z = Tot_z/length(Section_coords(:,1));
    
    Middle_point = [Av_x,Av_z];
end

% Function to reshape coordinate information for later use
function Coordinate_cell = Reshaper(Number_of_sections, Cavity_coordinates)
    Coordinate_cell = cell(Number_of_sections,1);
    Cavity_y_coordinates = sort(Cavity_coordinates(:,2));
    Cavity_y_coordinates_unique = unique(Cavity_y_coordinates);
    for a=1:Number_of_sections % Create separate section matrices for coordinates
        Section_index = 1; % Note: this is section i-1 as that starts counting at 0
        for b=1:length(Cavity_coordinates(:,1)) % Loop through the coordinate list
            if (Cavity_coordinates(b,2)==Cavity_y_coordinates_unique(a)) % Check y coordinates
                Vertex_matrix(a,Section_index,1) = Cavity_coordinates(b,1); % Store x coordinate
                Vertex_matrix(a,Section_index,2) = Cavity_coordinates(b,3); % Store z coordinate
                Section_index = Section_index + 1;
            end
        end
        clear b

        Section_coords = squeeze(Vertex_matrix(a,:,:)); % Extract the coordinates of the current section from the general matrix (we need 'squeeze' for that)
        for b=1:length(Section_coords(:,1))
            for c=1:length(Section_coords(1,:))
                if Section_coords(b,c) ~= 0 % The matrix has a lot of zeroes, since it must be filled according to its dimensions, and not all sections have the same number of coordinates
                    % Thus, thus check is needed to see how many vertices a
                    % section 'really' has
                    True_length = b;
                    True_width = c; % Sanity check; this must always be 2 (namely, x- and z-coordinate)
                end
            end
        end
        clear b c
        Section_coords = nonzeros(Section_coords); 
        Section_coords = reshape(Section_coords,True_length,True_width); % Remove the zeros from the section coordinates
        % Discussed earlier; what if an x- or z-coordinate is actually 0? Silly
        % trick: add 0.001 to it, so it doesn't get deleted

        % A way to store the section coordinates, but not the preferred way
        % NameLabelSection = genvarname(['Section_coords_',num2str(a)]);
        % eval([NameLabelSection '= Section_coords;']);
        % Why not preferred? It makes it more difficult to loop through these
        % coordinates

        % Better way: store them in a 'cell'
        Coordinate_cell{a} = Section_coords;
        clear Section_coords
    end
end

% Function to calculate the Euclidean distance between a middle point and
% the most outlying cavity point in the same section
function Max_distance = Distance_calculator(Middle_point,Coordinate_matrix)
    Max_distance = 0;
    for a=1:length(Coordinate_matrix(:,1))
        if sqrt((Coordinate_matrix(a,1)-Middle_point(1))^2+(Coordinate_matrix(a,2)-Middle_point(2))^2) > Max_distance
            Max_distance = sqrt((Coordinate_matrix(a,1)-Middle_point(1))^2+(Coordinate_matrix(a,2)-Middle_point(2))^2);
        end
    end
    clear a
end

% These voxels are always present, if we assume a block/beam model
function [Voxel,index] = Add_corner_voxels(index,Voxel,Limit_vector,Rock_density)
    
    % Bottom corner southwest
    Voxel(index,1) = Limit_vector(1) + 0.5;
    Voxel(index,2) = Limit_vector(3) + 0.5;
    Voxel(index,3) = Limit_vector(5) + 0.5;
    Voxel(index,4) = Rock_density;
    index = index + 1;
    
    % Bottom corner southeast
    Voxel(index,1) = Limit_vector(2) - 0.5;
    Voxel(index,2) = Limit_vector(3) + 0.5;
    Voxel(index,3) = Limit_vector(5) + 0.5;
    Voxel(index,4) = Rock_density;
    index = index + 1;
    
    % Bottom corner northwest
    Voxel(index,1) = Limit_vector(1) + 0.5;
    Voxel(index,2) = Limit_vector(4) - 0.5;
    Voxel(index,3) = Limit_vector(5) + 0.5;
    Voxel(index,4) = Rock_density;
    index = index + 1;
    
    % Bottom corner northeast
    Voxel(index,1) = Limit_vector(2) - 0.5;
    Voxel(index,2) = Limit_vector(4) - 0.5;
    Voxel(index,3) = Limit_vector(5) + 0.5;
    Voxel(index,4) = Rock_density;
    index = index + 1;
    
    % Top corner southwest
    Voxel(index,1) = Limit_vector(1) + 0.5;
    Voxel(index,2) = Limit_vector(3) + 0.5;
    Voxel(index,3) = Limit_vector(6) - 0.5;;
    Voxel(index,4) = Rock_density;
    index = index + 1;
    
    % Top corner southeast
    Voxel(index,1) = Limit_vector(2) - 0.5;
    Voxel(index,2) = Limit_vector(3) + 0.5;
    Voxel(index,3) = Limit_vector(6) - 0.5;
    Voxel(index,4) = Rock_density;
    index = index + 1;
    
    % Top corner northwest
    Voxel(index,1) = Limit_vector(1) + 0.5;
    Voxel(index,2) = Limit_vector(4) - 0.5;
    Voxel(index,3) = Limit_vector(6) - 0.5;
    Voxel(index,4) = Rock_density;
    index = index + 1;
    
    % Top corner northeast
    Voxel(index,1) = Limit_vector(2) - 0.5;
    Voxel(index,2) = Limit_vector(4) - 0.5;
    Voxel(index,3) = Limit_vector(6) - 0.5;
    Voxel(index,4) = Rock_density;
    index = index + 1;
    
end

% For the 'smart-grid' method, we initialise a grid consisting of only
% rock
function Starting_grid = Make_rocky_grid(Limit_vector,Rock_density,Grid_size,Grid)
    Matrix_size = Grid(1)*Grid(2)*Grid(3);
    
    Starting_grid = zeros(Matrix_size, 4);
    Element_index = 1;
    for x=1:Grid(1)
        for y=1:Grid(2)
            for z=1:Grid(3)
                Starting_grid(Element_index,1) = Limit_vector(1) + (x-1)*Grid_size;
                Starting_grid(Element_index,2) = Limit_vector(3) + (y-1)*Grid_size;
                Starting_grid(Element_index,3) = Limit_vector(5) + (z-1)*Grid_size;
                Starting_grid(Element_index,4) = Rock_density;
                
                Element_index = Element_index + 1;
            end
        end
    end
end

% A method to round the cavity coordinates to the existing grid
function Reduced_cavity_coordinates = Grid_cavity_coordinates(Cavity_coordinates,Grid_size,Starting_grid)
    for a=1:length(Cavity_coordinates(:,1))
        for b=1:length(Cavity_coordinates(1,:))
            [mincoor,index] = min(abs(Starting_grid(:,b)));
            Cavity_coordinates(a,b) = Starting_grid(index,b)+round(Cavity_coordinates(a,b)/Grid_size)*Grid_size;
        end
    end
    Reduced_cavity_coordinates = Cavity_coordinates;
end

% A function for the 'smart grid' method that calculates the 'dimension' of
% the cavity in the current section
function Max_distance = Grid_distance_calculator(Coordinate_array)
    Max_distance = 0;
    if length(Coordinate_array) == 1
        Max_distance = 0;
    else
        for a=1:length(Coordinate_array)
            for b=2:length(Coordinate_array)
                if abs(Coordinate_array(a)-Coordinate_array(b)) > Max_distance
                    Max_distance = abs(Coordinate_array(a)-Coordinate_array(b));
                end
            end
            clear b
        end
        clear a
    end
end

%% Function to convert a xml file to a MatLab structure
% From: Wouter Falkena
% (https://nl.mathworks.com/matlabcentral/fileexchange/28518-xml2struct)

function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

    if (nargin < 1)
        clc;
        help xml2struct
        return
    end
    
    if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
        % input is a java xml object
        xDoc = file;
    else
        %check for existance
        if (exist(file,'file') == 0)
            %Perhaps the xml extension was omitted from the file name. Add the
            %extension and try again.
            if (isempty(strfind(file,'.xml')))
                file = [file '.xml'];
            end
            
            if (exist(file,'file') == 0)
                error(['The file ' file ' could not be found']);
            end
        end
        %read the xml file
        xDoc = xmlread(file);
    end
    
    %parse xDoc into a MATLAB structure
    s = parseChildNodes(xDoc);
    
end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
    % Recurse over node children.
    children = struct;
    ptext = struct; textflag = 'Text';
    if hasChildNodes(theNode)
        childNodes = getChildNodes(theNode);
        numChildNodes = getLength(childNodes);

        for count = 1:numChildNodes
            theChild = item(childNodes,count-1);
            [text,name,attr,childs,textflag] = getNodeData(theChild);
            
            if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
                %XML allows the same elements to be defined multiple times,
                %put each in a different cell
                if (isfield(children,name))
                    if (~iscell(children.(name)))
                        %put existsing element into cell format
                        children.(name) = {children.(name)};
                    end
                    index = length(children.(name))+1;
                    %add new element
                    children.(name){index} = childs;
                    if(~isempty(fieldnames(text)))
                        children.(name){index} = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name){index}.('Attributes') = attr; 
                    end
                else
                    %add previously unknown (new) element to the structure
                    children.(name) = childs;
                    if(~isempty(text) && ~isempty(fieldnames(text)))
                        children.(name) = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name).('Attributes') = attr; 
                    end
                end
            else
                ptextflag = 'Text';
                if (strcmp(name, '#cdata_dash_section'))
                    ptextflag = 'CDATA';
                elseif (strcmp(name, '#comment'))
                    ptextflag = 'Comment';
                end
                
                %this is the text in an element (i.e., the parentNode) 
                if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                    if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                        ptext.(ptextflag) = text.(textflag);
                    else
                        %what to do when element data is as follows:
                        %<element>Text <!--Comment--> More text</element>
                        
                        %put the text in different cells:
                        % if (~iscell(ptext)) ptext = {ptext}; end
                        % ptext{length(ptext)+1} = text;
                        
                        %just append the text
                        ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                    end
                end
            end
            
        end
    end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
    % Create structure of node info.
    
    %make sure name is allowed as structure name
    name = toCharArray(getNodeName(theNode))';
    name = strrep(name, '-', '_dash_');
    name = strrep(name, ':', '_colon_');
    name = strrep(name, '.', '_dot_');

    attr = parseAttributes(theNode);
    if (isempty(fieldnames(attr))) 
        attr = []; 
    end
    
    %parse child nodes
    [childs,text,textflag] = parseChildNodes(theNode);
    
    if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
        %get the data of any childless nodes
        % faster than if any(strcmp(methods(theNode), 'getData'))
        % no need to try-catch (?)
        % faster than text = char(getData(theNode));
        text.(textflag) = toCharArray(getTextContent(theNode))';
    end
    
end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
    % Create attributes structure.

    attributes = struct;
    if hasAttributes(theNode)
       theAttributes = getAttributes(theNode);
       numAttributes = getLength(theAttributes);

       for count = 1:numAttributes
            %attrib = item(theAttributes,count-1);
            %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
            %attributes.(attr_name) = char(getValue(attrib));

            %Suggestion of Adrian Wanner
            str = toCharArray(toString(item(theAttributes,count-1)))';
            k = strfind(str,'='); 
            attr_name = str(1:(k(1)-1));
            attr_name = strrep(attr_name, '-', '_dash_');
            attr_name = strrep(attr_name, ':', '_colon_');
            attr_name = strrep(attr_name, '.', '_dot_');
            attributes.(attr_name) = str((k(1)+2):(end-1));
       end
    end
end