%% Function to calculate the maximum dimension of the current cavity section
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