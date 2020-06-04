
% Function to write a density to the voxels; either Rock_density or 0
function [Density,Count_index] = Calculate_density(Rock_density,Section_coordinates,Current_coordinate)

            % One of many ways to check whether a point is inside a 3D polyhedron.
            % Sample code obtained from https://nl.mathworks.com/matlabcentral/answers/101396-is-there-a-function-in-matlab-for-detecting-points-inside-a-polyhedron
            % 'Delaunayn' and 'Tsearchn' are MatLab standard funcions
            % Triangularization = delaunayn([Section_coordinates(:,1) Section_coordinates(:,2) Section_coordinates(:,3)]); % Generate delaunay triangulization
            % Search = tsearchn([Section_coordinates(:,1) Section_coordinates(:,2) Section_coordinates(:,3)], Triangularization, Current_coordinate); % Determine which triangle point is within
            % Inside_bool = ~isnan(Search); % Convert to logical vector
            
            % However, if the space to search through is 2D, the faster
            % approach 'inpolygon' can be used
            [In,On] = inpolygon(Current_coordinate(1),Current_coordinate(3),Section_coordinates(:,1),Section_coordinates(:,2));
            if (In==1 || On==1)
                Inside_bool = true;
            else
                Inside_bool = false;
            end
            
            if Inside_bool == true
                Density = 0;
            else
                Density = Rock_density;
            end

    % Maybe interesting alternative:
    % https://nl.mathworks.com/matlabcentral/answers/152189-volume-of-3d-polyhedron
end