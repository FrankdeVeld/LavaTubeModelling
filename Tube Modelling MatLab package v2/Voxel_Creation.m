% Function to transform a cavity model into a voxel model, filled with
% cavity voxels and rock voxels
function Voxel_Creation(Rock_density,Cavity_coordinates,Grid_matrix_3n,Cavity_bool,Model_limits,Grid_spacing,File_name)
    if Cavity_bool == true
        % The standard MatLab functionality of making 'alphaShapes' is used
        % here. If needed, they can be easily visualised with
        % 'plot(Cavity_Alphashape)'
        Cavity_Alphashape = alphaShape(Cavity_coordinates);
        % The reason for the use of alphashapes (earlier we already had the
        % full model) is the built-in efficient way of checking whether
        % points are in the interior. This is done for all points together
        % and is much faster than alternatives
        In_bool = inShape(Cavity_Alphashape,Grid_matrix_3n(:,1),Grid_matrix_3n(:,2),Grid_matrix_3n(:,3));
    else
        % Small trick: if there is no cavity to be modelled, just set every
        % point to 'false', i.e. not interior 
        In_bool = false(length(Grid_matrix_3n(:,1)),1);
    end
    
    % A way to sort the inside and outside points and to add the density as
    % a fourth column
    Inside_points = [Grid_matrix_3n(In_bool,1),Grid_matrix_3n(In_bool,2),Grid_matrix_3n(In_bool,3)];
    Inside_density = zeros(length(Inside_points(:,1)),1);

    Inside_points = [Inside_points,Inside_density];

    Outside_points = [Grid_matrix_3n(~In_bool,1),Grid_matrix_3n(~In_bool,2),Grid_matrix_3n(~In_bool,3)];
    Outside_density = Rock_density*ones(length(Outside_points(:,1)),1);

    Outside_points = [Outside_points,Outside_density];
    
    % 'Export' is the voxelised model and consists of interior and exterior
    % points with their respective densities.
    Export = [Inside_points;Outside_points];
    
    % This sorting is to have both files sorted the same way. It is useful
    % for checking, but should not affect performance
    [~,idx] = sort(Export(:,3)); % sort just the first column
    Export = Export(idx,:); 
    [~,idx] = sort(Export(:,2)); % sort just the first column
    Export = Export(idx,:); 
    [~,idx] = sort(Export(:,1)); % sort just the first column
    Export = Export(idx,:); 
   
    % It was suggested to, instead of adding edging voxels, just add corner
    % voxels and count on IGMAS interpolation. This is expected to be
    % speeding up the process, but so far it has not been proven to work;
    % thus it is commented out
    
%   Export = Add_corner_voxels(Export, Rock_density, Model_limits, Grid_spacing);
    
    % The actual process of writing the voxel file
    fileID = fopen(File_name,'w');
    fprintf(fileID,'x\t y\t z\t cellValue\n');
    for i = 1:length(Export(:,1))
      fprintf(fileID,'%d\t%d\t%d\t%d\n', Export(i,1),Export(i,2),Export(i,3),Export(i,4));
    end
    fclose(fileID);
end

% Function to add corner voxels in a beam-shaped model. The function itself
% works, but the interpolation in IGMAS does not behave as expected.In the
% current script, this function is not used
function Export = Add_corner_voxels(Export, Rock_density, Model_limits, Grid_spacing)
    Base_length = length(Export(:,1));

    for a=1:8
        if mod(a,2) == 1
            Export(Base_length+a,1) = Model_limits(2)-Grid_spacing/2;
            Export(Base_length+a,2) = Model_limits(3)+Grid_spacing/2;
        else
            Export(Base_length+a,1) = Model_limits(1)+Grid_spacing/2;
            Export(Base_length+a,2) = Model_limits(4)-Grid_spacing/2;
        end
        if a<5
            Export(Base_length+a,3) = Model_limits(5)+Grid_spacing/2;
        else
            Export(Base_length+a,3) = Model_limits(6)-Grid_spacing/2;
        end
        Export(Base_length+a,4) = Rock_density;
    end
end