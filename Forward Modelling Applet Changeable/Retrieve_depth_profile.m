%% Goal of this function: provide the dpeth profile of the cavity to the applet
% Author: Frank de Veld
% Date: June 22th, 2020
% Output: A matrix of depth, position coordinates

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
    Coordinates = Coordinates - max(Coordinates(:,3));
    Outline = boundary(Coordinates(:,1),Coordinates(:,2));
    
    x_coords = Coordinates(:,1);
    y_coords = Coordinates(:,2);
    z_coords = Coordinates(:,3);
    Outline_coordinates = [x_coords(Outline),y_coords(Outline),z_coords(Outline)];
    
    [Sorted_coords_x, Sorted_coords_y] = poly2cw(Outline_coordinates(:,1), Outline_coordinates(:,2));
    
    New_sorted_coords_x = [];
    New_sorted_coords_y = [];
    
    if (strcmp(Lava_tube,'Name1')==true)
        for b=1:length(Outline_coordinates(:,1))
            New_index = b+430; % Empirical to make it start at one of the ends
            New_index = mod(New_index,length(Outline_coordinates(:,1)))+1;
            New_sorted_coords_x(1,b) = Sorted_coords_x(New_index);
            New_sorted_coords_y(1,b) = Sorted_coords_y(New_index);
        end
    elseif (strcmp(Lava_tube,'Name2')==true)
        for b=1:length(Outline_coordinates(:,1))
            New_index = b+624; % Empirical to make it start at one of the ends
            New_index = mod(New_index,length(Outline_coordinates(:,1)))+1;
            New_sorted_coords_x(1,b) = Sorted_coords_x(New_index);
            New_sorted_coords_y(1,b) = Sorted_coords_y(New_index);
        end
    elseif (strcmp(Lava_tube,'Name3')==true)
        for b=1:length(Outline_coordinates(:,1))
            New_index = b+48; % Empirical to make it start at one of the ends
            New_index = mod(New_index,length(Outline_coordinates(:,1)))+1;
            New_sorted_coords_x(1,b) = Sorted_coords_x(New_index);
            New_sorted_coords_y(1,b) = Sorted_coords_y(New_index);
        end
    end
    
    Halfway_point = floor(length(Outline_coordinates(:,1))/2);
    
    z_array = [];
    Distance_array = [];
    Distance_array(1) = 0;
    Middle_points_matrix = [];
    figure()
    hold on
    for b=1:Halfway_point
        Point_1 = [New_sorted_coords_x(b),New_sorted_coords_y(b)];
        Point_2 = [New_sorted_coords_x(length(Outline_coordinates(:,1))-b),New_sorted_coords_y(length(Outline_coordinates(:,1))-b)];
        plot([Point_1(1),Point_2(1)],[Point_1(2),Point_2(2)])
        
        [~, Index_1]=ismember(Point_1,[Outline_coordinates(:,1),Outline_coordinates(:,2)],'rows');
        z_1 = Outline_coordinates(Index_1,3);
        
        [~, Index_2]=ismember(Point_2,[Outline_coordinates(:,1),Outline_coordinates(:,2)],'rows');
        z_2 = Outline_coordinates(Index_2,3);
        
        z_array(b) = (z_1+z_2)/2;
        
        Middle_points_matrix(b,1) = (Point_1(1) + Point_2(1))/2;
        Middle_points_matrix(b,2) = (Point_1(2) + Point_2(2))/2;
        
        if b>1
            Distance_array(b) = sqrt((Middle_points_matrix(b,1)-Middle_points_matrix(b-1,1))^2+(Middle_points_matrix(b,2)-Middle_points_matrix(b-1,2))^2) + Distance_array(b-1);
        end
    end

    figure()
    plot(Distance_array,z_array)
    
    Output = [Distance_array;z_array]';
    
    Depth_profile_name = ['Depth_profile_',Lava_tube];
    save(Depth_profile_name,'Output');
end


figure()
plot(alphaShape(Coordinates))