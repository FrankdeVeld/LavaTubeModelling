%% Function to write .vxo-files
% Densities are filled in according to the known cavity location
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

    Number_samples_per_section = 5;
    Starting_grid = Make_rocky_grid(Limit_vector+Grid_size/2,Rock_density,Grid_size,Grid);
    
    if Cavity_bool == true
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

                    Temp_density = Calculate_density(Rock_density,Current_section_matrix,Temp_coord); % Density

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
    end
    Voxel = Starting_grid;


    % Writing it to a .vxo file
    fileID = fopen(File_name,'w');
    fprintf(fileID,'x\t y\t z\t cellValue\n');
    for i = 1:length(Voxel(:,1))
      fprintf(fileID,'%d\t%d\t%d\t%d\n', Voxel(i,1),Voxel(i,2),Voxel(i,3),Voxel(i,4));
    end
    fclose(fileID);
end


%% Auxilary Functions

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