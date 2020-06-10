%% This script is meant for making models in MatLab for IGMAS
% Authors: Frank de Veld, Stephanie Bringeland
% Date: May 27th, 2020
% Output: .model file

% Goal: take a set of coordinates of a cave as input and
% have a .model file as output without side effect issues

% Advantages: easy to change, fast both in producing the files and importing
% them in IGMAS
% Disadvantages: limited resolution, script is difficult to understand, per
% section, only one cavity can be reliably modelled.

% Later task: find the best way of sectioning
% Later task: Fix the annoying ³-bug
% Later task: Make it so that multiple cavities can be modelled
% Later task: look at the error one makes when approximating by comparing
% the true coordinates with the processed coordinates

clear all
close all

% What can be changed?
% - Directory of the lava tube model file (Step 1)
% - Reducing factor: how many data points do you find acceptable (Step 2)
% - Vertices per section: how many points do you want per section (Step 2)
% - Alternatively: how many sections do you want (Step 2)
% For these, note that at least three points are needed in every section.
% Taylor the numbers to the number of unique x or y-coordinates in the
% cavity model
% - Cavity depth: The distance between cavity and stations (Step 3)
% - Rock density: what the density of the surrounding rock is (Step 3)
% - Station resolution: how many stations do you want on the grid (Step 3)
% - Name of the model (Step 3)

% Later task: also make the size of the cavity changeable 

%% Step 1: Import cavity coordinates

% Either: import a file directly displaying coordinates 
Cavity_coordinates = importdata('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Matlab model maker\Cavity_Test_Data_Coords.csv'); % Change the directory accordingly

% Or: import a 3D model file (for example a lava tube model), and take one
% of its properties to be the coordinate file
% [A,~] = stlread('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Lava tube models\cueva-de-los-siete-lagos\source\agua_text\Los lagos ok.stl');
% Cavity_coordinates = A.Points;
%% Step 2: Do some data processing on these coordinates

% First: reduce the number of points (for example, by a factor 50)
% Idea: each 50th coordinate is picked for a fair representation of points
Reduce_factor = 1;

Reduced_index = 1;
Reduced_cavity_coordinates = zeros(floor(length(Cavity_coordinates(:,1))/Reduce_factor),3);
for a = 1:length(Cavity_coordinates(:,1))
    if mod(a,Reduce_factor)==0
        Reduced_cavity_coordinates(Reduced_index,:) = Cavity_coordinates(a,:);
        Reduced_index = Reduced_index + 1;
    end
end
clear a

% The 'Cavity_coordinates' variable is reset
Cavity_coordinates = Reduced_cavity_coordinates*100;

% Count occurrence of all y-coordinates, on which we do the sectioning
y_coords = Cavity_coordinates(:,2);

% A fairly inefficient way of checking whether there are enough vertices
% per section (enough being 3 or more). However, this has not been a major
% issue so far
Bin_count = [1:length(Cavity_coordinates(:,1))];
y_occur = zeros(length(y_coords),1);
index = 1;
for a=1:length(y_coords)
    for b=1:length(y_coords) % Nested loops; very brute-force...
        if y_coords(a) == y_coords(b)
            y_occur(index) = y_occur(index) + 1; 
        end
    end
    index = index + 1;
end
clear a b

Interpolation_bool = false; % Indicates whether interpolation is needed
for a=1:length(y_occur)
    if y_occur(a) < 3 % Less than three vertices per section will not work
        Interpolation_bool = true;
    end
end

% Perhaps a preferred way is to just say the following for those large
% files
% Interpolation_bool = true

% Later task: don't restructure the whole data set if only a few
% y-coordinates have occurrences less than 3.

if Interpolation_bool == true
    % If this is the case, a restructure of the data set must take place
    % before it can be processed. This whole piece of code restrucures the
    % data such that at least 3 vertices are present per section.
    
    % Two ways of interpolation (more can be made)
    Fixed_number_vertices_bool = true; % Method 1: have a fixed number of vertices per section
    Fixed_number_section_bool = false; % Method 2: have a fixed number of sections
    
    if Fixed_number_vertices_bool==true % Method 1: fixed number of vertices per section
        % Choose how many points you want per section:
        Vertices_per_section = 3;

        Number_of_sections = floor(length(Cavity_coordinates(:,1))/Vertices_per_section); % Flooring, as you want at a minimum three vertices per section
        Leftovers = mod(length(Cavity_coordinates(:,1)),Vertices_per_section); % For if the section number is not a divisor of the number of coordinates
        Leftovers_front = ceil(Leftovers/2); % Crude method: put half of the left-overs in the first section
        Leftovers_back = floor(Leftovers/2); % Crude method: put half of the left-overs in the last section

        [~,Sorted_y] = sort(Cavity_coordinates(:,2)); % Sort the y-coordinates
        Sorted_cavity_coordinates = Cavity_coordinates(Sorted_y,:); % Sort the whole matrix according to these y-coordinates

        for a=1:Leftovers_front % The 'leftovers' are partly placed in the first section
            Middle_coor = Leftovers_front + ceil(Vertices_per_section/2); % This is the y-coordinate of the middle of a group of coordinates of size 'vertices_per_section'
            Sorted_cavity_coordinates(a,2) = Sorted_cavity_coordinates(Middle_coor,2); % Let these coordinates 'adopt' the one of the middle of the group
        end
        clear a

        Current_section_index = Leftovers_front + 1; % The first few coordinates already have been fixed; we start at 'Leftovers_front + 1'
        for a=1:Number_of_sections % Loop over the sections you want (defined earlier)
            for b=1:Vertices_per_section % For each section, we want this number of vertices
                Middle_coor = Leftovers_front + (a-1)*Vertices_per_section + ceil(Vertices_per_section/2); % Previously to this section, the 'leftovers_front' and (a-1) sections have been gone through
                Sorted_cavity_coordinates(Current_section_index,2) = Sorted_cavity_coordinates(Middle_coor,2); % Let all coordinates of the section adopt the same y-coordinate 
                Current_section_index = Current_section_index + 1; % So that the indexing of 'Sorted_cavity_coordinates' still works
            end
        end
        clear a

        for a=(length(Sorted_cavity_coordinates(:,2))-Leftovers_back):length(Sorted_cavity_coordinates(:,2)) % The same thing for the back happens as what has been done to the front. Note the start and end
            Middle_coor = length(Sorted_cavity_coordinates(:,2)) - Leftovers_back - floor(Vertices_per_section/2);
            Sorted_cavity_coordinates(a,2) = Sorted_cavity_coordinates(Middle_coor,2);
        end
        clear a

        Cavity_coordinates = Sorted_cavity_coordinates; % I am not sure whether this is needed, but it's useful for checking
    elseif Fixed_number_section_bool == true % Method 2: fixed number of sections in total
        
        Number_of_sections = 100; % It is advised to taylor this number to the total number of coordinates.
        
        y_range = max(Cavity_coordinates(:,2))-min(Cavity_coordinates(:,2)); % The range is needed to fairly distribute the sections.
        
        Leftover = mod(length(Cavity_coordinates(:,2)),(length(Cavity_coordinates(:,2))/Number_of_sections)); % Idea: assign x points to one section, but with this you have leftovers
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
            y_coor = min(Cavity_coordinates(:,2)) + 0.5*y_range/Number_of_sections; % Thus, the leftovers at the front are combined in the first section
            Sorted_cavity_coordinates(a,2) = y_coor;
        end
        clear a
        
        for a=(Leftover_front+1):(length(Sorted_cavity_coordinates(:,2))-Leftover_back-1) % This is a loop through the sorted coordinates
            Current_section = floor(a/(length(Cavity_coordinates(:,1))/Number_of_sections)); % The information on the section we are working in is still needed
            y_coor = min(Cavity_coordinates(:,2)) + y_range/Number_of_sections*Current_section + 0.5*y_range/Number_of_sections;
            % Explanation: this is the y-value, so the minimum is
            % min(Cavity_coordinates(:,2)) and the middle of the current
            % section is half a section away from the 'nodes' explained
            % before. Lastly, y_range/Number_of_sections*Current_section is
            % to account for previous sections
            Sorted_cavity_coordinates(a,2) = y_coor;
        end
        clear a
        
        for a=(length(Sorted_cavity_coordinates(:,2))-Leftover_back):length(Sorted_cavity_coordinates(:,2)) % The same thing happens for the back as happened for the front
            y_coor = min(Cavity_coordinates(:,2)) + y_range - 0.5*y_range/Number_of_sections;
            Sorted_cavity_coordinates(a,2) = y_coor;
        end
        clear a
        
        Cavity_coordinates = Sorted_cavity_coordinates;
    end
end

% A bit of a silly limitation; zero values for x or z-coordinates are not
% 'allowed' as they give problems when removing zeros from matrices. Since
% quite large approximations are already present, it doesn't matter so much
% if another 0.001 is added.
for a=1:length(Cavity_coordinates(:,1))
    if Cavity_coordinates(a,1) == 0
        Cavity_coordinates(a,1) = Cavity_coordinates(a,1) + 0.001;
    end
    if Cavity_coordinates(a,3) == 0
        Cavity_coordinates(a,3) = Cavity_coordinates(a,3) + 0.001;
    end
end
% If this is a major problem later on, this change can be reversed later on
% in the file

% An idea I had; round the coordinates for more interpolation. It doesn't
% turn out to be needed
% Cavity_coordinates = round(Cavity_coordinates,3);

%% Step 3: define parameters, do set-up
% Define important parameters
Cavity_depth = 1; % depth in meter
Rock_density = 2.5; % density in g/cm³
Sizing_factor = 0.1; % To deal with edge effects
Model_name = ['Test_Large']; % Model_name = ['Los_Lagos_May27_reduce10_100sections']; % The name you want your model to have
Stations_bool = true; % if you want to also export stations
Station_resolution = 50; % if you want stations, how many do you want (this is the square root of it)

%% Step 4: open up the file, define the bodies
feature('DefaultCharacterSet', 'UTF8'); % Needed mostly to make the ³-character work...
File_name = ['Sec_',Model_name,'.model']; % The file extension is always .model
fid = fopen(File_name,'w'); % Open the file (this also works if the file does not exist yet)
fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n'); % Make sure the encoding of the file works
fprintf(fid, '<geodata name="Test Model">\n'); % Open the geodata class; the main class of the whole file

% Define a reference body (never needs to be changed, apart from maybe the colour)
fprintf(fid, '<property name="body" value="reference"> \n    <property name="density" units="g/cm³" value="0.0"></property> \n <color red="0.5019608" green="0.5019608" blue="0.5019608"></color> \n</property> \n' );

% Define another body. Note: density can be changed!
String_rock = ['<property name="body" value="stone"> \n    <property name="density" units="g/cm³" value="%4.2f"></property> \n <color red="1" green="0" blue="0"></color> \n</property> \n'];
String_rock = sprintf(String_rock,Rock_density); % This is so that a generic name can be put in the line above, and the value can be put here
fprintf(fid, String_rock);

% Define a cavity. Assumed density: 0
fprintf(fid, '<property name="body" value="cavity"> \n    <property name="density" units="g/cm³" value="0"></property> \n <color red="0" green="0" blue="1"></color> \n</property> \n' );

% Define 'wider' model to circumvent issues with boundary effects, based on
% the dimensions of the cavity
% Here, we need to look how big these extensions must be
% Interpretation: the total model has the size of the cavity, times 2*Sizing_factor.
% This can be tweaked!
Limit_vector = [];

Cavity_limit_vector = [];
for a=1:3
    Cavity_limit_vector = [Cavity_limit_vector,min(Cavity_coordinates(:,a)),max(Cavity_coordinates(:,a))];
end

Cavity_length_x = Cavity_limit_vector(2)-Cavity_limit_vector(1);
Cavity_length_y = Cavity_limit_vector(4)-Cavity_limit_vector(3);

Limit_vector(1) = round(Cavity_limit_vector(1) - Sizing_factor*Cavity_length_x);
Limit_vector(2) = round(Cavity_limit_vector(2) + Sizing_factor*Cavity_length_x);
Limit_vector(3) = round(Cavity_limit_vector(3) - Sizing_factor*Cavity_length_y);
Limit_vector(4) = round(Cavity_limit_vector(4) + Sizing_factor*Cavity_length_y);
Limit_vector(5) = round(Cavity_limit_vector(5) - Cavity_depth);
Limit_vector(6) = round(Cavity_limit_vector(6) + Cavity_depth);

Cavity_coordinates(:,3) = Cavity_coordinates(:,3) - Limit_vector(6);
Limit_vector(5)  = Limit_vector(5) - Limit_vector(6);
Limit_vector(6)  = Limit_vector(6) - Limit_vector(6);


%% Step 5: export a .model file without cavity

% Export a 'blank' model for comparison
Export_no_cavity_model(Rock_density, Model_name, Limit_vector);

%% Step 6: create sections, restructure the cavity data for later use
% For now: sections parallel to x axis, thus perpendicular to y-axis
% Later task: find out what the optimal way of sectioning is

Cavity_y_coordinates = sort(Cavity_coordinates(:,2));
Cavity_y_coordinates_unique = unique(Cavity_y_coordinates); % Look for all unique y coordinates. At this point, each unique y-coordinate should have at least two neighbouring points with equal y-coordinates
Number_of_sections = length(Cavity_y_coordinates_unique); % For each unique y-coordinate, a section ought to be present

Vertex_cell = cell(Number_of_sections,1);
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
    Vertex_cell{a} = Section_coords;
    clear Section_coords
end
clear a 

%% Step 7: fill in the sections with vertices

% The .model file starts with a 0-vertex and a 0-section
Vertex_number = 0;
Section_number = 0;
Visualisation_factor = 1;

% Create first section for preventing edge effects
Section_coords = [];
[Vertex_number, Section_number] = Section_writing(Section_coords,Limit_vector,Vertex_number,Section_number,fid,Cavity_y_coordinates_unique,true,Visualisation_factor);

% Create the cavity sections
for a=1:Number_of_sections
    Section_coords = cell2mat(Vertex_cell(a,1)); % Cell2mat is needed to read the coordinates correctly
    [Vertex_number, Section_number] = Section_writing(Section_coords,Limit_vector,Vertex_number,Section_number,fid,Cavity_y_coordinates_unique,false,Visualisation_factor);
end
clear a

% Create last section for preventing edge effects
[Vertex_number, Section_number] = Section_writing(Section_coords,Limit_vector,Vertex_number,Section_number,fid,Cavity_y_coordinates_unique,true,Visualisation_factor);

%% Step 8: close the file

fprintf(fid, '</geodata>\n' );  
fclose(fid);  % Close the geodata class, and with this the file

%% Step 9: if needed, also export stations
if (Stations_bool==true)
    Export_stations_file(Limit_vector,Station_resolution,Model_name) % ASSUMING flat land, also currently only a set grid is created without other patterns
end

%% Step 10: if needed, also export a model where the cavity is not inside the sections. Result: nice visualisation of the cavity
% This step was a result of a sign error somewhere, placing the cavity
% above z=0 instead of underneath. However, this error actually provides us
% with one of the best ways to look at the cavity in 3D and how it has been
% modelled, which is nice
Visualisation_bool = true;
if Visualisation_bool == true
    Visualisation_factor = -1;
    Export_3D_model_file(Limit_vector,Rock_density,Number_of_sections,Vertex_cell,Cavity_y_coordinates_unique,Model_name,Visualisation_factor) 
end

%% Step 11: import the models and station file in IGMAS, calculate the anomalies, export the results

%(The goal is to automate this)

%% Step 12: import the results from IGMAS and make a grid of the difference in observed signal (e.g. gZ)

% Import the .stations file with results. The 'xml2struct' function is a
% function from file exchange
% (https://nl.mathworks.com/matlabcentral/fileexchange/28518-xml2struct)
% which makes it easy to extract data (coordinates, values) from files like
% the .stations file.

% Put this at true if you ahve results to analyse
Calculation_bool = false;

if Calculation_bool == true

    With_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\TestsMay27\Sec_With_cavity_results_Los_Lagos_May27_reduce10_100sections.stations');

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

    No_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\TestsMay27\Sec_No_cavity_results_Los_Lagos_May27_reduce10_100sections.stations');

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
    figure()
    scatter3(With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3))
    title('Difference between the signals of models with and without a cavity','interpreter','latex')
    xlabel('x-coordinate [m]','interpreter','latex')
    ylabel('y-coordinate [m]','interpreter','latex')
    zlabel('Gravity difference [mGal]','interpreter','latex')

end

%% Functions for specific tasks

% Function to write the sections and their vertices in an output file
function [Vertex_number, Section_number] = Section_writing(Section_coords,Limit_vector,Vertex_number,Section_number,fid,Unique_y_coordinates,Edge_sections_boolean,Visualisation_factor)
    % The 'Visualisation_factor' is only checked for the optional
    % visualisation part, where just a '3D-model' of the cave is exported,
    % without the surrounding rock

    % 'Temp_text' is a string that has all the information of the current
    % section, and at the end of this function this text gets added to the
    % .model-file
    Temp_text = [];
    if (Edge_sections_boolean==true) % Edge sections for preventing edge effects. These never have cavities
        if Vertex_number == 0
            Which_side = 3;
        else
            Which_side = 4;
        end
        
        Temp_text = ['<geometry>' ...
          '<cross_section name="%d" '...
          'x_start="%4.2f" y_start="%4.2f" x_end="%4.2f" y_end="%4.2f"></cross_section> \n'];
        Temp_text = sprintf(Temp_text,Section_number,Limit_vector(1),Limit_vector(Which_side),Limit_vector(2),Limit_vector(Which_side)); % The x-limits are defined by the body extent. The y-limit is the minimum y-coordinate

        Vertex_number_old = Vertex_number; % An 'old' variable is defined to not lose track of the starting value in this section 

        Base_vertex = 1;
        for a=1:4 % Base Vertices. We assume four vertices, placed at the edged of the extent of the body, as defined earlier
            [String_appendix, Vertex_number] = Vertex_writing(Section_coords,Limit_vector,Vertex_number,false,Base_vertex,Visualisation_factor);
            Temp_text = [Temp_text,String_appendix];
            Base_vertex = Base_vertex+1; % The order matters; it is clockwise, and it is useful to keep track of this
        end
        clear a
        
        Temp_text = [Temp_text,'<entry type="polygon" id_list="%d %d %d %d "> \n' ...
             '<property name="body" value="stone"></property> \n'...
          '</entry> \n']; % Definition of the base body; the stone. Always consisting of four points, in this version
        Temp_text = sprintf(Temp_text,Vertex_number_old,Vertex_number_old+1,Vertex_number_old+2,Vertex_number_old+3);

        Temp_text = [Temp_text,'</geometry> \n'];
        fprintf(fid,Temp_text);
    else
        Vertex_number_old = Vertex_number;

            Temp_text = ['<geometry>' ...
              '<cross_section name="%d" '...
              'x_start="%4.2f" y_start="%4.2f" x_end="%4.2f" y_end="%4.2f"></cross_section> \n'];
            Temp_text = sprintf(Temp_text,Section_number,Limit_vector(1),Unique_y_coordinates(Section_number),Limit_vector(2),Unique_y_coordinates(Section_number));

        if Visualisation_factor == 1   
            Base_vertex = 1;
            for a=1:4 % Base Vertices
                [String_appendix, Vertex_number] = Vertex_writing(Section_coords,Limit_vector,Vertex_number,false,Base_vertex,Visualisation_factor);
                Temp_text = [Temp_text,String_appendix];
                Base_vertex = Base_vertex+1;
            end
            clear a
        end

        Middle_point = Middle_point_seeker(Section_coords); % Finds the middle point of all vertices. Output: [x,z]-coordinate
        Sorted_coords = Clockwise_vertex_sorter(Section_coords,Middle_point); % Sorts the vertices such that they are clockwise. Output: (n,2)-matrix

        Cavity_vertex = 1;
        for a=1:length(Section_coords(:,1)) % Cavity Vertices
            [String_appendix, Vertex_number] = Vertex_writing(Sorted_coords,Limit_vector,Vertex_number,true,Cavity_vertex,Visualisation_factor);
            Temp_text = [Temp_text,String_appendix];
            Cavity_vertex = Cavity_vertex+1;
        end
        clear a
        
        % Add the base body definition to the file
        if Visualisation_factor == 1
            Temp_text = [Temp_text,'<entry type="polygon" id_list="%d %d %d %d "> \n' ...
                 '<property name="body" value="stone"></property> \n'...
              '</entry> \n']; % Definition of the base body; the stone. Always consisting of four points, in this version
            Temp_text = sprintf(Temp_text,Vertex_number_old,Vertex_number_old+1,Vertex_number_old+2,Vertex_number_old+3);

            Vertex_number_after_base = Vertex_number_old + 3; % To keep track; after the base vertices have been added, this is the starting vertex number
        end
        if Visualisation_factor == -1
            Vertex_number_after_base = Vertex_number_old-1;
        end
        % Definition of the cavity. Since it can consist of an arvitrary number
        % of points, it is a bit trickier to define, but it does the same
        % as the base-vertex writing previously. Note that since it is
        % unknown in advance how many nodes there are, writing the id's
        % must be done in a loop.
        Temp_text = [Temp_text,'<entry type="polygon" id_list="']; 
        for a=1:length(Sorted_coords(:,1)) 
            Temp_text = [Temp_text,'%d '];
            Temp_text = sprintf(Temp_text,Vertex_number_after_base + a);
        end
        clear a
        
        % Finish the definition of the cavity body
        Temp_text = [Temp_text,'"> \n', ...
             '<property name="body" value="cavity"></property> \n'...
          '</entry> \n'];

        % Finish the writing of this section, and add it to the file
        Temp_text = [Temp_text,'</geometry> \n'];
        fprintf(fid,Temp_text);
    end
    Section_number = Section_number+1;
end

% This function is called when making a section, and addresses the wiritng
% of vertices. Note that these vertices must be created clockwise in a
% section
function [Section_text,Vertex_number] = Vertex_writing(Section_coords,Limit_vector,Vertex_number,Cavity_vertex,Order_vertex,Visualisation_factor)   
    % This way for the base vertices always works, it seems (working with a
    % block, thas is). Here also, it is assumed that the base can be
    % defined by four vertices in a rectangle
    if Visualisation_factor == 1
        if (Cavity_vertex == false) % Base vertices
            Section_text = ['<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'];
            if Order_vertex == 1
                Section_text = sprintf(Section_text,Vertex_number,Limit_vector(2)-Limit_vector(1),Limit_vector(5)-Limit_vector(6));
            elseif Order_vertex == 2
                Section_text = sprintf(Section_text,Vertex_number,0,Limit_vector(5)-Limit_vector(6));
            elseif Order_vertex == 3
                Section_text = sprintf(Section_text,Vertex_number,0,0);
            elseif Order_vertex == 4
                Section_text = sprintf(Section_text,Vertex_number,Limit_vector(2)-Limit_vector(1),0);   
            end
        end
    end
    
    % Note: coordinates are with respect to the limits of the section
    % limits, and the second coordinate is the distance to the first
    % coordinate, not to the 'origin'.
    if (Cavity_vertex == true) % Cavity vertices. The order has already been defined by the clock-wise function earlier.
        Section_text = ['<vertex id="%d" x="%4.2f" z="%4.2f"></vertex> \n'];
        Section_text = sprintf(Section_text,Vertex_number,Section_coords(Order_vertex,1)-Limit_vector(1),Section_coords(Order_vertex,2)); % Filling in the vertex number and coordinates
    end
    
    Vertex_number = Vertex_number + 1;
end

% Function to find the middle of a cavity in a section; useful for ordering
% vertices
function Middle_point = Middle_point_seeker(Section_coords) 
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

% Function to sort points such that the are ordered clockwise; this is needed
% as IGMAS reads the points in clockwise order.
function Sorted_coords = Clockwise_vertex_sorter(Section_coords,Middle_point) 
    Angle_conversion = [];
    Complex_conversion = [];
    for a=1:length(Section_coords(:,1))
        Complex_conversion(a) = (Section_coords(a,1)-Middle_point(1)) + 1i*(Section_coords(a,2)-Middle_point(2)); % The points are converted to complex numbers with as origin the middle point
        Angle_conversion(a) = angle(Complex_conversion(a)); % With this conversion, the complex argument gives a way to order the points
    end
    clear a
    Sorted_angles = sort(Angle_conversion,'descend'); % Clockwise means descending arguments
    for a=1:length(Section_coords(:,1))
        for b=1:length(Section_coords(:,1))
            if Angle_conversion(a) == Sorted_angles(b)
                Section_coords_new(b,1) = Section_coords(a,1); % Create a sorted list
                Section_coords_new(b,2) = Section_coords(a,2);
            end
        end
    end
    clear a b
    Sorted_coords = Section_coords_new;
end

% Function to also export a model without cavity. Since it is assumed that
% this base is just a beam, it is fairly easy to define and the code is
% just a repetition of earlier code
function Export_no_cavity_model(Rock_density, Model_name, Limit_vector)
    feature('DefaultCharacterSet', 'UTF8'); % Needed mostly to make the ³-character work...
    Model_name = ['Sec_No_cavity_',Model_name,'.model'];
    fid = fopen(Model_name,'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n'); % I am not sure what this does; result in proper encoding?
    fprintf(fid, '<geodata name="Test Model">\n'); % Open the geodata class

    % Define a reference body
    fprintf(fid, '<property name="body" value="reference"> \n    <property name="density" units="g/cm³" value="0.0"></property> \n <color red="0.5019608" green="0.5019608" blue="0.5019608"></color> \n</property> \n' );

    % Define another body. Note: density can be changed!
    String_rock = ['<property name="body" value="stone"> \n    <property name="density" units="g/cm³" value="%4.2f"></property> \n <color red="1" green="0" blue="0"></color> \n</property> \n'];
    String_rock = sprintf(String_rock,Rock_density);
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
    Temp_text = sprintf(Temp_text,Limit_vector(2)-Limit_vector(1),Limit_vector(5)-Limit_vector(6),0,Limit_vector(5)-Limit_vector(6),0,0,Limit_vector(2)-Limit_vector(1),0);
    Section_text = [Section_text,Temp_text];

    Section_text = [Section_text,'<entry type="polygon" id_list="0 1 2 3 "> \n' ...
         '<property name="body" value="stone"></property> \n'...
      '</entry> \n']; 

    Section_text = [Section_text,'</geometry> \n'];
    fprintf(fid,Section_text);
    
    % Section 1
    
    Section_text = ['<geometry>' ...
      '<cross_section name="1" '...
      'x_start="%4.2f" y_start="%4.2f" x_end="%4.2f" y_end="%4.2f"></cross_section> \n'];
    Section_text = sprintf(Section_text,Limit_vector(1),Limit_vector(4),Limit_vector(2),Limit_vector(4));

    Temp_text = ['<vertex id="4" x="%4.2f" z="%4.2f"></vertex> \n' ...
    '<vertex id="5" x="%4.2f" z="%4.2f"></vertex> \n' ....
    '<vertex id="6" x="%4.2f" z="%4.2f"></vertex> \n'...
    '<vertex id="7" x="%4.2f" z="%4.2f"></vertex> \n'];
    Temp_text = sprintf(Temp_text,Limit_vector(2)-Limit_vector(1),Limit_vector(5),0,Limit_vector(5),0,Limit_vector(6),Limit_vector(2)-Limit_vector(1),Limit_vector(6));
    Section_text = [Section_text,Temp_text];

    Section_text = [Section_text,'<entry type="polygon" id_list="4 5 6 7 "> \n' ...
         '<property name="body" value="stone"></property> \n'...
      '</entry> \n']; 

    Section_text = [Section_text,'</geometry> \n'];
    fprintf(fid,Section_text);


    fprintf(fid, '</geodata>\n' );  % Close the geodata class, and with this the file
    fclose(fid);
end

% Function to also export a stations file. There is no elevation assumed,
% and also just a standard grid. This can be improved in the future.
function Export_stations_file(Limit_vector,Station_resolution, Model_name)
    feature('DefaultCharacterSet', 'UTF8'); % Needed mostly to make the ³-character work...
    File_name = ['Sec_Stations_',Model_name,'.stations'];
    fid = fopen(File_name,'w');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n'); % I am not sure what this does; result in proper encoding?
    Temp_text = ['<geodata> \n' ...
    '<projection name="unknown" units="m"></projection> \n']; % Open the geodata class
    fprintf(fid, Temp_text);
    x_increment = (Limit_vector(2)-Limit_vector(1))/Station_resolution;
    y_increment = (Limit_vector(4)-Limit_vector(3))/Station_resolution;
    for a=1:(Station_resolution+1)
        for b=1:(Station_resolution+1)
            Temp_text = ['<vertex x="%4.2f" y="%4.2f" z="%4.2f"> \n</vertex> \n'];
            Temp_text = sprintf(Temp_text,Limit_vector(1) + x_increment*(a-1),Limit_vector(3) + y_increment*(b-1),Limit_vector(6));
            fprintf(fid,Temp_text);
        end
    end
    clear a b
    fprintf(fid, '</geodata>');
end

function Export_3D_model_file(Limit_vector,Rock_density,Number_of_sections,Vertex_cell,Cavity_y_coordinates_unique,Model_name,Visualisation_factor)
    
    feature('DefaultCharacterSet', 'UTF8'); % Needed mostly to make the ³-character work...
    File_name = ['Sec_Visualisation_',Model_name,'.model'];
    fid = fopen(File_name,'w'); % Open the file (this also works if the file does not exist yet)
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n'); % I am not sure what this does; result in proper encoding?
    fprintf(fid, '<geodata name="Test Model">\n'); % Open the geodata class; the main class of the whole file

    % Define a reference body (never needs to be changed, apart from maybe the colour)
    fprintf(fid, '<property name="body" value="reference"> \n    <property name="density" units="g/cm³" value="0"></property> \n <color red="0.5019608" green="0.5019608" blue="0.5019608"></color> \n</property> \n' );

    % Define another body. Note: density can be changed!
    String_rock = ['<property name="body" value="stone"> \n    <property name="density" units="g/cm³" value="%4.2f"></property> \n <color red="1" green="0" blue="0"></color> \n</property> \n'];
    String_rock = sprintf(String_rock,Rock_density); % This is so that a generic name can be put in the line above, and the value can be put here
    fprintf(fid, String_rock);

    % Define a cavity. Assumed density: 0
    fprintf(fid, '<property name="body" value="cavity"> \n    <property name="density" units="g/cm³" value="0"></property> \n <color red="0" green="0" blue="1"></color> \n</property> \n' );


    % The .model file starts with a 0-vertex and a 0-section
    Vertex_number = 0;
    Section_number = 1;

    
    % Create the cavity sections
    for a=1:Number_of_sections
        Section_coords = cell2mat(Vertex_cell(a,1)); % Cell2mat is needed to read the coordinates correctly
        [Vertex_number, Section_number] = Section_writing(Section_coords,Limit_vector,Vertex_number,Section_number,fid,Cavity_y_coordinates_unique,false,Visualisation_factor);
    end
    clear a

    fprintf(fid, '</geodata>\n' );  
    fclose(fid);  % Close the geodata class, and with this the file
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