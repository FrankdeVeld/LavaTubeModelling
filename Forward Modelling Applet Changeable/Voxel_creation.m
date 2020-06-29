%% Goal of this function: transform a set of coordinates signifying the boundaries of a cavity to a voxel script with the cavity filled
% Author: Frank de Veld
% Date: June 17th, 2020
% Output:  one .vxo file

% Function to transform a cavity model into a voxel model, filled with
% cavity voxels and rock voxels

% This function is called automatically when running the 'Voxel_main.m'
% script
function Voxel_creation(Model_density,Cavity_coordinates,Grid_matrix_3n,File_name)
    % The standard MatLab functionality of making 'alphaShapes' is used
    % here. If needed, they can be easily visualised with
    % 'plot(Cavity_Alphashape)'
    Cavity_Alphashape = alphaShape(Cavity_coordinates);
    Cavity_Alphashape.Alpha = Cavity_Alphashape.Alpha*2; % Some security which might have to be improved later on
    % The reason for the use of alphashapes (earlier we already had the
    % full model) is the built-in efficient way of checking whether
    % points are in the interior. This is done for all points together
    % and is much faster than alternatives
    In_bool = inShape(Cavity_Alphashape,Grid_matrix_3n(:,1),Grid_matrix_3n(:,2),Grid_matrix_3n(:,3));
    
    % A way to sort the inside and outside points and to add the density as
    % a fourth column
    Inside_points = [Grid_matrix_3n(In_bool,1),Grid_matrix_3n(In_bool,2),Grid_matrix_3n(In_bool,3)];
    Inside_density = Model_density*ones(length(Inside_points(:,1)),1);

    Inside_points = [Inside_points,Inside_density];

    Outside_points = [Grid_matrix_3n(~In_bool,1),Grid_matrix_3n(~In_bool,2),Grid_matrix_3n(~In_bool,3)];
    Outside_density = zeros(length(Outside_points(:,1)),1);

    Outside_points = [Outside_points,Outside_density];
    
    % 'Export' is the voxelised model and consists of interior and exterior
    % points with their respective densities.
    Export = [Inside_points;Outside_points];
    
    if(length(Inside_points(:,1)) < 1)
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

    % The actual process of writing the voxel file
    fileID = fopen(File_name,'w');
    fprintf(fileID,'x\t y\t z\t cellValue\n');
    for i = 1:length(Export(:,1))
      fprintf(fileID,'%d\t%d\t%d\t%d\n', Export(i,1),Export(i,2),Export(i,3),Export(i,4));
    end
    fclose(fileID);
end