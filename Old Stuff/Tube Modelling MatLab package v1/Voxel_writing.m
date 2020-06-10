%% Function to write .vxo-files
% Densities are filled in according to the known cavity location
function [Vox_limits] =  Voxel_writing(Rock_density,Cavity_coordinates,Cavity_limit_vector,Cavity_bool,File_name,Grid_spacing,Number_voxels_per_direction,Vox_limits,Sizing_factor,Model_name)
    Voxel = []; % Initialize the full data file
   
    

    Starting_grid = Make_rocky_grid(Cavity_limit_vector,Rock_density,Grid_spacing,Number_voxels_per_direction);
    
    if Cavity_bool == true
       
        Gridded_cavity_coordinates = Grid_cavity_coordinates(Cavity_coordinates,Grid_spacing,Starting_grid);
    
        Gridded_cavity_limit_vector = [];
        for a=1:3
            Gridded_cavity_limit_vector = [Gridded_cavity_limit_vector,round(min(Gridded_cavity_coordinates(:,a))),round(max(Gridded_cavity_coordinates(:,a)))];
        end
        
        Unique_gridded_cavity_coordinates = unique(Gridded_cavity_coordinates, 'rows');

        [~,idx] = sort(Unique_gridded_cavity_coordinates(:,2)); % sort just the first column
        Sorted_unique_gridded_cavity_coordinates = Unique_gridded_cavity_coordinates(idx,:);   % sort the whole matrix using the sort indices
        
        
        
        xrange = Gridded_cavity_limit_vector(2) - Gridded_cavity_limit_vector(1);
        yrange = Gridded_cavity_limit_vector(4) - Gridded_cavity_limit_vector(3);
        zrange = Gridded_cavity_limit_vector(6) - Gridded_cavity_limit_vector(5);

        Gridded_cavity_range_vector = [xrange,yrange,zrange];

        Vox_min_x = min(Gridded_cavity_coordinates(:,1)) - Gridded_cavity_range_vector(1)*Sizing_factor;
        Vox_max_x = max(Gridded_cavity_coordinates(:,1)) + Gridded_cavity_range_vector(1)*Sizing_factor;
        Vox_min_y = min(Gridded_cavity_coordinates(:,2)) - Gridded_cavity_range_vector(2)*Sizing_factor;
        Vox_max_y = max(Gridded_cavity_coordinates(:,2)) + Gridded_cavity_range_vector(2)*Sizing_factor;
        Vox_min_z = min(Gridded_cavity_coordinates(:,3)) - Gridded_cavity_range_vector(3)*Sizing_factor;
        Vox_max_z = 0;
        
        Vox_min_x = Vox_min_x + mod((Vox_max_x-Vox_min_x),Grid_spacing);
        Vox_min_y = Vox_min_y + mod((Vox_max_y-Vox_min_y),Grid_spacing);
        Vox_min_z = Vox_min_z + mod((Vox_max_z-Vox_min_z),Grid_spacing);

        Vox_limits = [Vox_min_x,Vox_max_x,Vox_min_y,Vox_max_y,Vox_min_z,Vox_max_z];
        
        % Note: this model has almost nothing; the proper cavity will be defined
        % later with voxels
        Write_basic_model_file(Rock_density, Model_name, Vox_limits, Cavity_limit_vector, Number_voxels_per_direction(2)-1); % -1, as the number of intermediate sections is required as input
        
        % Due to technical issues, we rather don't have zeros in the
        % matrix. If this is really a problem, a method can be implemented
        % that temporarily changes the zero coordinates to another value.
        for a=1:length(Sorted_unique_gridded_cavity_coordinates(:,1))
            if Sorted_unique_gridded_cavity_coordinates(a,1) == 0
                Sorted_unique_gridded_cavity_coordinates(a,1) = Sorted_unique_gridded_cavity_coordinates(a,1) + 0.001;
            end
            if Sorted_unique_gridded_cavity_coordinates(a,2) == 0
                Sorted_unique_gridded_cavity_coordinates(a,2) = Sorted_unique_gridded_cavity_coordinates(a,2) + 0.001;
            end
            if Sorted_unique_gridded_cavity_coordinates(a,3) == 0
                Sorted_unique_gridded_cavity_coordinates(a,3) = Sorted_unique_gridded_cavity_coordinates(a,3) + 0.001;
            end
        end

        Sorted_gridded_cavity_y_coordinates = Sorted_unique_gridded_cavity_coordinates(:,2);
        Sorted_unique_gridded_cavity_y_coordinates = unique(Sorted_gridded_cavity_y_coordinates);
        Number_cavity_sections = length(Sorted_unique_gridded_cavity_y_coordinates);
        Coordinate_cell = Reshaper(Number_cavity_sections,Sorted_unique_gridded_cavity_coordinates);

        index = 1;
        
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

            Number_x_point_evaluation = Max_x_distance/Grid_spacing + 1;
            Number_z_point_evaluation = Max_z_distance/Grid_spacing + 1;

            Min_x = min(Current_section_matrix(:,1));
            Min_z = min(Current_section_matrix(:,2));

            for b=1:Number_x_point_evaluation
                for c=1:Number_z_point_evaluation
                    Temp_x = Min_x + b*Grid_spacing;
                    Temp_z = Min_z + c*Grid_spacing;
                    Temp_y = Sorted_unique_gridded_cavity_y_coordinates(a);

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
    
    Voxel = Add_corner_voxels(Voxel, Rock_density, Grid_spacing, Vox_limits);

    % Writing it to a .vxo file
    fileID = fopen(File_name,'w');
    fprintf(fileID,'x\t y\t z\t cellValue\n');
    for i = 1:length(Voxel(:,1))
      fprintf(fileID,'%d\t%d\t%d\t%d\n', Voxel(i,1),Voxel(i,2),Voxel(i,3),Voxel(i,4));
    end
    fclose(fileID);
end


%% Auxilary Functions

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

function [Voxel,Vox_limits] = Add_corner_voxels(Voxel, Rock_density, Grid_spacing, Vox_limits)
    Base_length = length(Voxel(:,1));
    
    
    for a=1:8
        if mod(a,2) == 1
            Voxel(Base_length+a,1) = Vox_limits(2)-Grid_spacing/2;
            Voxel(Base_length+a,2) = Vox_limits(3)+Grid_spacing/2;
        else
            Voxel(Base_length+a,1) = Vox_limits(1)+Grid_spacing/2;
            Voxel(Base_length+a,2) = Vox_limits(4)-Grid_spacing/2;
        end
        if a<5
            Voxel(Base_length+a,3) = Vox_limits(5)+Grid_spacing/2;
        else
            Voxel(Base_length+a,3) = Vox_limits(6)-Grid_spacing/2;
        end
        Voxel(Base_length+a,4) = Rock_density;
    end
end