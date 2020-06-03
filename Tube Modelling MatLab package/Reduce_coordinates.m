%% Function to reduce the number of coordinates in a file in the most even way
% The function takes every xth coordinate and saves that, x being the
% 'Reduce_factor'

function Reduced_cavity_coordinates = Reduce_coordinates(Cavity_coordinates, Reduce_factor)
    Reduced_index = 1;
    Reduced_cavity_coordinates = zeros(floor(length(Cavity_coordinates(:,1))/Reduce_factor),3);
    for a = 1:length(Cavity_coordinates(:,1))
        if mod(a,Reduce_factor)==0
            Reduced_cavity_coordinates(Reduced_index,:) = Cavity_coordinates(a,:);
            Reduced_index = Reduced_index + 1;
        end
    end
    clear a
end