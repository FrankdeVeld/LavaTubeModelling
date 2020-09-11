%% Goal of this script: make a video of improvement when certain measurements are done
% Authors: Frank de Veld
% Date: September 8th, 2020
% Output: something visual

Model_name = 'LAGOS';
Model_multiplication = 10;
Model_density = 2.25; 
Cavity_depth = 40; 
Number_of_points = 500;
Number_of_revs = 10;
% This comes from the voxel main file

% Systematically creating the files
MM = Model_multiplication;
MD = Model_density;
CD = Cavity_depth;

Station_spacing = 2.5*MM;

Model_limits = [-1255,3281,-1038,2709];

Limits_x = 4536;
Limits_y = 3747;

Centre(1) = Model_limits(1) + Limits_x/2;
Centre(2) = Model_limits(3) + Limits_y/2;

if Limits_x > Limits_y
    End_point(1) = Centre(1);
    End_point(2) = Model_limits(4);
else
    End_point(1) = Model_limits(2);
    End_point(2) = Centre(2);
end

Number_of_stations(1) = 182;
Number_of_stations(2) = 151;




% For 'xml2struct', see https://nl.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
File_to_open = ['.\Database Result Files\Results\',Model_name,'_MM_',num2str(MM),'_MD_',num2str(MD),'_CD_',num2str(CD),'_results.stations'];
Results = xml2struct(File_to_open);

Results_array = [];
Z = [];
xindex = 1;
yindex = 1;
Station_number_x = Number_of_stations(1); 
Station_number_y = Number_of_stations(2); 
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

    Z(xindex,yindex) = str2num(Temp_result);
    yindex = yindex + 1;
end
x = unique(Results_array(:,1));
y = unique(Results_array(:,2));
[X,Y] = meshgrid(x,y);
X = X';
Y = Y';

Z_spiral = NaN*ones(length(X(:,1)),length(X(1,:)));

NoP = Number_of_points;
NoR = Number_of_revs;

x_values = X(:,1);
Positive_indices_x = x_values > 0;
min_x = min(x_values(Positive_indices_x));

y_values = Y(1,:);
Positive_indices_y = y_values > 0;
min_y = min(y_values(Positive_indices_y));

Results_array_shift = [Results_array(:,1) - min_x, Results_array(:,2) - min_y];
%            Results_array_only_xy = [Results_array(:,1),Results_array(:,2)];
Spiral_points = Gridded_spiral(Centre, End_point, NoP, NoR, Station_spacing)';

Spiral_points_ordered(1:25,:) = flipud(Spiral_points(26:50,:));
Spiral_points_ordered(26:50,:) = flipud(Spiral_points(1:25,:));

%%
MovieMatrix = struct('cdata',[],'colormap',[]);
Z_spiral = NaN*ones(length(X(:,1)),length(X(1,:)));

v = VideoWriter('10Turns500PointsWithCavity.avi');
open(v);

File_name = ['.\Data Analysis\Outline_LAGOS.mat'];
Outline = MM*importdata(File_name);

for d=1:length(Spiral_points_ordered(:,1))
    [tf, index]=ismember(Spiral_points_ordered(d,:),Results_array_shift,'rows');

    Spiral_index_x = ceil(index/Number_of_stations(2));

    if mod(index,Station_number_y+1)==0
        Spiral_index_y = Station_number_y;
    else
        Spiral_index_y = mod(index, Station_number_y+1) ;
    end

    Z_spiral(Spiral_index_x,Spiral_index_y) = Z(Spiral_index_x,Spiral_index_y);

    Z_spiral_int = inpaint_nans(Z_spiral,3);

    if d==1
       base = Z_spiral_int; 
    end

    if mod(d,1) == 0
%         figure()
%         contour(X,Y,Z_spiral_int-base)
        
        figure()
        hold on
        
        Contour_levels = 0:0.001:2;
        h1 = contour(X,Y,Z_spiral_int,Contour_levels);
        h2 = scatter(Spiral_points_ordered(1:d,1),Spiral_points_ordered(1:d,2),25,'k','o','filled');
        h3 = plot(Outline(:,1),Outline(:,2),'r--');
        caxis([0,2])
        c = colorbar;
        c.Label.String = 'Gravity anomaly (mGal)';
        
        str = sprintf('Contour map after %d gravity measurements, 10 turns',d);
        title(str,'interpreter','latex')
        xlabel('x-coordinate [m]','interpreter','latex')
        ylabel('y-coordinate [m]','interpreter','latex')
        legend([h2,h3],'Measurement points','Cavity outline', 'Location', 'northwest')
        hold off
        
        frame = getframe(gcf);
        for t=1:2
            writeVideo(v,frame);
        end
        
        MovieMatrix(:,round(d/1)) = getframe(gcf);
%         figure()
%         contour(X,Y,Z_spiral)
    end
end   

close(v)

File_name_spiral = ['.\Database Result Files\Results\',Model_name,'_MM_',num2str(MM),'_MD_',num2str(MD),'_CD_',num2str(CD),'_NR_',num2str(NoR),'_NP_',num2str(NoP),'_spiral.mat'];
if isfile(File_name_spiral)==false
    save(File_name_spiral,'Z_spiral')
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