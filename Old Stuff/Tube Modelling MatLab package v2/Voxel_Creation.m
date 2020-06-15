% Function to transform a cavity model into a voxel model, filled with
% cavity voxels and rock voxels
function Voxel_Creation(Rock_density,Cavity_density,Cavity_coordinates,Grid_matrix_3n,Cavity_bool,Model_limits,Grid_spacing,File_name,Visualisation_bool)
    if Cavity_bool == true
        % The standard MatLab functionality of making 'alphaShapes' is used
        % here. If needed, they can be easily visualised with
        % 'plot(Cavity_Alphashape)'
        Cavity_Alphashape = alphaShape(Cavity_coordinates);
        Cavity_Alphashape.Alpha = Cavity_Alphashape.Alpha*2; % Some security which might has to be improved later on
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
    Inside_density = Cavity_density*ones(length(Inside_points(:,1)),1);

    Inside_points = [Inside_points,Inside_density];
    
    if Visualisation_bool == true
        Vis_name = ['Visualisation_',File_name];
        fileID = fopen(Vis_name,'w');
        fprintf(fileID,'x\t y\t z\t cellValue\n');
        for i = 1:length(Inside_points(:,1))
          fprintf(fileID,'%d\t%d\t%d\t%d\n', Inside_points(i,1),Inside_points(i,2),Inside_points(i,3),Inside_points(i,4));
        end
        fclose(fileID);
    end
    
    Outside_points = [Grid_matrix_3n(~In_bool,1),Grid_matrix_3n(~In_bool,2),Grid_matrix_3n(~In_bool,3)];
    Outside_density = Rock_density*ones(length(Outside_points(:,1)),1);

    Outside_points = [Outside_points,Outside_density];
    
    % 'Export' is the voxelised model and consists of interior and exterior
    % points with their respective densities.
    Export = [Inside_points;Outside_points];
    
    if(Cavity_bool == true && length(Inside_points(:,1)) < 1)
        error('Error: no inside points detected. Decrease grid spacing, or increase alpha of alphaShape')
    end
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
    
  Export = Add_corner_voxels(Export, Rock_density, Model_limits, Grid_spacing);
    
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
        if (a == 1 || a == 5)
            Export(Base_length+a,1) = Model_limits(1)+Grid_spacing/2;
            Export(Base_length+a,2) = Model_limits(3)+Grid_spacing/2;
        elseif (a==2 || a == 6)
            Export(Base_length+a,1) = Model_limits(1)+Grid_spacing/2;
            Export(Base_length+a,2) = Model_limits(4)-Grid_spacing/2;
        elseif (a==3 || a == 7)
            Export(Base_length+a,1) = Model_limits(3)-Grid_spacing/2;
            Export(Base_length+a,2) = Model_limits(2)+Grid_spacing/2;
        elseif (a==4 || a == 8)
            Export(Base_length+a,1) = Model_limits(3)-Grid_spacing/2;
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