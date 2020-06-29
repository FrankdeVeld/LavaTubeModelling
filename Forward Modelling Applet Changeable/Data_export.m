%% Goal of this script: transform an IGMAS results file into X, Y and Z matrices for later quick visualisation
% Authors: Frank de Veld
% Date: June 17th, 2020
% Output:  one to three .mat files with results and station coordinates

% Run when you have results of all files (of a certain lava tube model), as
% this script will make results in a loop. Takes about 4 seconds per
% configuration, so certainly when you have about 100 configurations per
% lava tube model you can take a quick break here

% This should be known beforehand
Lava_tube_model_array = {'Name1','Name2','Name3','Name4'};
Model_multiplication_array = [1, 2, 5, 10]; % Leads to cavities of about 20 m diameter to 200 m diameter
Model_density_array = [1.5, 2.25, 3]; % Densities are easy to interpret and interpolate
Cavity_depth_array = [10, 20, 40, 80]; % Such that most cavities are still structurally stable
Station_spacing_array = [1, 2, 3, 4, 5, 10, 50]; 
Number_of_points_array = [25, 50, 100, 150, 250];
Number_of_revs_array = [2,3,4,5,6];
% This comes from the voxel main file

File_directory_LAGOS = '.\Model files\Name1';
Cave_info = stlread(File_directory_LAGOS);
Cavity_coordinates_clean = Cave_info.vertices;
Edging_factor = 0.50;
Lava_tube_model = Lava_tube_model_array{1};

Station_spacing_vector = [];

% Systematically creating the files
Lava_tube_model = Lava_tube_model_array{1};
for a=3:4 % Sizing
    for b=1:3 % Model Density 
        for c=1:4 % Cavity Depth
            MM = Model_multiplication_array(a);
            MD = Model_density_array(b);
            CD = Cavity_depth_array(c);
            
            Station_spacing = 2.5*MM;
            Grid_spacing = 2.5*MM;
            
            Cavity_coordinates = MM*Cavity_coordinates_clean;
            
            Sizing_addition = [];
            Model_limits = [];
            for d=1:3 % In this loop (3 or 3x2 coordinates) the grid is created and the limits for the cavity and the model (with edging factor) are calculated
                Gridded_coords(d) = {round(min(Cavity_coordinates(:,d))):Grid_spacing:(round(max(Cavity_coordinates(:,d)))-Grid_spacing)};
                Sizing_addition(d) = round((Gridded_coords{1,d}(end)-Gridded_coords{1,d}(1))*Edging_factor/Grid_spacing)*Grid_spacing;
                if d==1
                    Model_limits(2*d-1) = round(min(Cavity_coordinates(:,d)))-Sizing_addition(d);
                    Model_limits(2*d) = round(max(Cavity_coordinates(:,d)))+Sizing_addition(d);
                elseif d ==2
                    Model_limits(2*d-1) = round(min(Cavity_coordinates(:,d)))-Sizing_addition(d); % Can be + for 'cut-off' cavity
                    Model_limits(2*d) = round(max(Cavity_coordinates(:,d)))+Sizing_addition(d); % Can be - for 'cut-off' cavity
                end
            end


            Limits_x = Model_limits(2) - Model_limits(1);
            Limits_y = Model_limits(4) - Model_limits(3);
            
            Centre(1) = Model_limits(1) + Limits_x/2;
            Centre(2) = Model_limits(3) + Limits_y/2;
            
            if Limits_x > Limits_y
                End_point(1) = Centre(1);
                End_point(2) = Model_limits(4);
            else
                End_point(1) = Model_limits(2);
                End_point(2) = Centre(2);
            end
            
            Station_spacing_vector(1) = round(Limits_x/Station_spacing)+1;
            Station_spacing_vector(2) = round(Limits_y/Station_spacing)+1;
            
            Results_writer(Lava_tube_model,MM,MD,CD,Station_spacing_vector,Station_spacing_array,Number_of_points_array,Number_of_revs_array,Centre,End_point,Station_spacing)
        end
    end
end

function Results_writer(Model_name,MM,MD,CD,Station_spacing_vector,Station_spacing_array,Number_of_points_array,Number_of_revs_array,Centre,End_point,Station_spacing)
    % For 'xml2struct', see https://nl.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
    File_to_open = ['.\Database generation\Results\',Model_name,'_MM_',num2str(MM),'_MD_',num2str(MD),'_CD_',num2str(CD),'_results.stations'];
    Results = xml2struct(File_to_open);

    Results_array = [];
    Z = [];
    xindex = 1;
    yindex = 1;
    Station_number_x = Station_spacing_vector(1); 
    Station_number_y = Station_spacing_vector(2); 
    for a=1:length(Results.geodata.vertex(1,:))
        % This trickery is needed to make the long list of coordinates into a
        % square matrix
        xindex = ceil(a/Station_number_y);
        if mod(yindex, Station_number_y)==0
            yindex = Station_number_y;
        else
            yindex = mod(yindex, Station_number_y) ;
        end

        Temp_x = Results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = Results.geodata.vertex{1,a}.Attributes.y;

        % Calculated gz
        Temp_result = Results.geodata.vertex{1,a}.property.Attributes.value;

        Results_array(a,1) = str2num(Temp_x);
        Results_array(a,2) = str2num(Temp_y);
        Results_array(a,3) = str2num(Temp_result);
        
        Saving_array_3D(a,1) = str2num(Temp_x);
        Saving_array_3D(a,2) = str2num(Temp_y);
        Saving_array_3D(a,3) = str2num(Temp_result);

        Z(xindex,yindex) = str2num(Temp_result);
        yindex = yindex + 1;
    end
    x = unique(Results_array(:,1));
    y = unique(Results_array(:,2));
    [X,Y] = meshgrid(x,y);
    X = X';
    Y = Y';
    
    Z_spiral = NaN*ones(length(X(:,1)),length(X(1,:)));
    
    for b=1:length(Number_of_points_array)
        for c=1:length(Number_of_revs_array)
            NoP = Number_of_points_array(b);
            NoR = Number_of_revs_array(c);
            
            x_values = X(:,1);
            Positive_indices_x = x_values > 0;
            min_x = min(x_values(Positive_indices_x));
            
            y_values = Y(1,:);
            Positive_indices_y = y_values > 0;
            min_y = min(y_values(Positive_indices_y));
            
            Results_array_shift = [Results_array(:,1) - min_x, Results_array(:,2) - min_y];
%            Results_array_only_xy = [Results_array(:,1),Results_array(:,2)];
            Spiral_points = Gridded_spiral(Centre, End_point, NoP, NoR, Station_spacing)';
            
            Results_array_spiral_3D = [];
            
            for d=1:length(Spiral_points(:,1))
                [tf, index]=ismember(Spiral_points(d,:),Results_array_shift,'rows');
                
                Spiral_index_x = ceil(index/Station_spacing_vector(2));
                
                if mod(index,Station_spacing_vector(2))==0
                    Spiral_index_y = Station_number_y;
                else
                    Spiral_index_y = mod(index, Station_number_y) ;
                end
                
                Z_spiral(Spiral_index_x,Spiral_index_y) = Z(Spiral_index_x,Spiral_index_y);
                
                Results_array_spiral_3D(d,:) = [Spiral_points(d,:),Z(Spiral_index_x,Spiral_index_y)];
            end   
            
            File_name_spiral = ['.\Database generation\Results\',Model_name,'_MM_',num2str(MM),'_MD_',num2str(MD),'_CD_',num2str(CD),'_NR_',num2str(NoR),'_NP_',num2str(NoP),'_spiral.mat'];
            if isfile(File_name_spiral)==false
                save(File_name_spiral,'Z_spiral')
            end
            
            File_name_spiral_3D = ['.\Database generation\Results\',Model_name,'_MM_',num2str(MM),'_MD_',num2str(MD),'_CD_',num2str(CD),'_NR_',num2str(NoR),'_NP_',num2str(NoP),'_spiral_3D.mat'];
            if isfile(File_name_spiral_3D)==false
                save(File_name_spiral_3D,'Results_array_spiral_3D')
            end
        end
    end
    
    
    % Reasoning: the X and Y coordinates of the stations are equal for all
    % configurations of the same lava tube model (apart from additional sizing). Thus, they need to be
    % saved only once per model
    for b=1:length(Station_spacing_array)
        X_SSF = X(1:Station_spacing_array(b):end,1:Station_spacing_array(b):end);
        Y_SSF = Y(1:Station_spacing_array(b):end,1:Station_spacing_array(b):end);
        Z_SSF = Z(1:Station_spacing_array(b):end,1:Station_spacing_array(b):end);
        
        Results_array_3D_SSF = [];
        Results_index = 1;
        for c=1:length(X_SSF(:,1))
            for d=1:length(X_SSF(1,:))
                Results_array_3D_SSF(Results_index,1) = X_SSF(c,d);
                Results_array_3D_SSF(Results_index,2) = Y_SSF(c,d);
                Results_array_3D_SSF(Results_index,3) = Z_SSF(c,d);
                
                Results_index = Results_index + 1;
            end
        end

        %Results_array_3D_SSF = Results_array(1:Station_spacing_array(b):end,:);
        
        File_name_X = ['.\Database generation\Results\',Model_name,'_MM_',num2str(MM),'_SSF_',num2str(Station_spacing_array(b)),'_X.mat'];
        if isfile(File_name_X)==false
            save(File_name_X,'X_SSF')
        end
        File_name_Y = ['.\Database generation\Results\',Model_name,'_MM_',num2str(MM),'_SSF_',num2str(Station_spacing_array(b)),'_Y.mat'];
        if isfile(File_name_Y)==false
            save(File_name_Y,'Y_SSF')
        end

        File_name_Z = ['.\Database generation\Results\',Model_name,'_MM_',num2str(MM),'_MD_',num2str(MD),'_CD_',num2str(CD),'_SSF_',num2str(Station_spacing_array(b)),'_Z.mat'];
        if isfile(File_name_Z)==false
            save(File_name_Z,'Z_SSF')
        end

        File_name_3D = ['.\Database generation\Results\',Model_name,'_MM_',num2str(MM),'_MD_',num2str(MD),'_CD_',num2str(CD),'_SSF_',num2str(Station_spacing_array(b)),'_3D.mat'];
        if isfile(File_name_3D)==false
            save(File_name_3D,'Results_array_3D_SSF');
        end
    end
end

function Spiral_points = Gridded_spiral(Centre, End_point, Number_of_points, Number_of_revs, Station_spacing)
    %random 2d points

    %radius of first point to second
    r = norm(Centre-End_point);
    %angle between two point wrt the y-axis
    if (Centre(1)-End_point(1))==0
        theta_offset=pi/2;
    else
        theta_offset = tan((End_point(2)- Centre(2))/(End_point(1)-Centre(1)));
    end

    t = linspace(0,r,Number_of_points); %radius as spiral decreases
    theta = linspace(0,2*pi*Number_of_revs,Number_of_points) + theta_offset; %angle as spiral decreases
    x = cos(theta).*t+Centre(1); 
    y = sin(theta).*t+Centre(2);

    x_gridded = round((x)./Station_spacing)*Station_spacing;
    y_gridded = round((y)./Station_spacing)*Station_spacing;

    Spiral_points = [x_gridded;y_gridded];
end