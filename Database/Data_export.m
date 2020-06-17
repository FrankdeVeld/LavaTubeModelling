%% Goal of this script: transform an IGMAS results file into X, Y and Z matrices for later quick visualisation
% Authors: Frank de Veld
% Date: June 17th, 2020
% Output:  one to three .mat files with results and station coordinates

% Run when you have results of all files (of a certain lava tube model), as
% this script will make results in a loop. Takes about 4 seconds per
% configuration, so certainly when you have about 100 configurations per
% lava tube model you can take a quick break here

% This should be known beforehand
Lava_tube_model_array = ['LAGOS','GALA','JAMEO','ETNA'];
Grid_spacing_array = [0.5, 1, 2.5, 5, 10];
Model_density_array = [1.5, 2, 2.5, 3, 3.5];
Cavity_depth_array = [1, 5, 10, 25, 50];

% This comes from the voxel main file
Station_resolution = 50;

% Systematically creating the files
Lava_tube_model = Lava_tube_model_array(1);
for a=1:5 % Grid Spacing
    for b=1:5 % Model Density 
        for c=1:5 % Cavity Depth
            GS = Grid_spacing_array(a);
            MD = Model_density_array(b);
            CD = Cavity_depth_array(c);
            Results_writer(Lava_tube_model,GS,MD,CD,Station_resolution)
        end
    end
end

function Results_writer(Model_name,GS,MD,CD,Station_resolution)
    % For 'xml2struct', see https://nl.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
    File_to_open = ['.\Database generation\Results\',Model_name,'_GS_',num2str(GS),'_MD_',num2str(MD),'_CD_',num2str(CD),'_results.stations'];
    Results = xml2struct(File_to_open);

    Results_array = [];
    Z = [];
    xindex = 1;
    yindex = 1;
    Station_number = Station_resolution+1; % Of course dependent on the station resolution, which usually is 50
    for a=1:length(Results.geodata.vertex(1,:))
        % This trickery is needed to make the long list of coordinates into a
        % square matrix
        xindex = ceil(a/Station_number);
        if mod(yindex, Station_number)==0
            yindex = Station_number;
        else
            yindex = mod(yindex, Station_number) ;
        end

        Temp_x = Results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = Results.geodata.vertex{1,a}.Attributes.y;

        % Calculated gz
        Temp_result = Results.geodata.vertex{1,a}.property.Attributes.value;

        Results_array(a,1) = str2num(Temp_x);
        Results_array(a,2) = str2num(Temp_y);
        Results_array(a,3) = str2num(Temp_result);

        Z(xindex,yindex) = str2num(Temp_result);
        yindex = yindex + 1;
    end
    clear a
    x = unique(Results_array(:,1));
    y = unique(Results_array(:,2));
    [X,Y] = meshgrid(x,y);
    X = X';
    Y = Y';

    % Reasoning: the X and Y coordinates of the stations are equal for all
    % configurations og the same lava tube model. Thus, they need to be
    % saved only once per model
    File_name_X = ['.\Database generation\Results\',Model_name,'_X.mat'];
    if isfile(File_name_X)==false
        save(X,File_name_X)
    end
    File_name_Y = ['.\Database generation\Results\',Model_name,'_Y.mat'];
    if isfile(File_name_Y)==false
        save(Y,File_name_Y)
    end
    
    File_name_Z = ['.\Database generation\Results\',Model_name,'_GS_',num2str(GS),'_MD_',num2str(MD),'_CD_',num2str(CD),'_Z.mat'];
    save(File_name_Z,'Z');
end