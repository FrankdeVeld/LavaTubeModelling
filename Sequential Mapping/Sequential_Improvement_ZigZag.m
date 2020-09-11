%% Goal of this script: make a video of improvement when certain measurements are done
% Authors: Frank de Veld
% Date: September 8th, 2020
% Output: something visual

Model_name = 'LAGOS';
Model_multiplication = 10;
Model_density = 2.25; 
Cavity_depth = 40; 
Number_of_points = 50;
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



x_values = X(:,1);
Positive_indices_x = x_values > 0;
min_x = min(x_values(Positive_indices_x));

y_values = Y(1,:);
Positive_indices_y = y_values > 0;
min_y = min(y_values(Positive_indices_y));

Results_array_shift = [Results_array(:,1) - min_x, Results_array(:,2) - min_y];
%            Results_array_only_xy = [Results_array(:,1),Results_array(:,2)];


x_grid = Model_limits(1)+250:round((Model_limits(2)-Model_limits(1))/7):Model_limits(2)-250;
y_grid = Model_limits(3)+250:round((Model_limits(4)-Model_limits(3))/7):Model_limits(4)-250;

[X_grid,Y_grid] = meshgrid(x_grid,y_grid);

X_grid = round((X_grid)./Station_spacing)*Station_spacing;
Y_grid = round((Y_grid)./Station_spacing)*Station_spacing;

Plot_array = [];

for a=1:length(x_grid)
   for b=1:length(y_grid)
        Plot_array(length(y_grid)*(a-1)+b,1) = X_grid(a,b);
        Plot_array(length(y_grid)*(a-1)+b,2) = Y_grid(a,b);
   end
end

Sorted_plot_array = [];
Sorted_plot_array(1:7,:) = Plot_array(1:7,:);
Sorted_plot_array(8:14,:) = flip(Plot_array(8:14,:));
Sorted_plot_array(15:21,:) = Plot_array(15:21,:);
Sorted_plot_array(22:28,:) = flip(Plot_array(22:28,:));
Sorted_plot_array(29:35,:) = Plot_array(29:35,:);
Sorted_plot_array(36:42,:) = flip(Plot_array(36:42,:));
Sorted_plot_array(43:49,:) = Plot_array(43:49,:);

%%
MovieMatrix = struct('cdata',[],'colormap',[]);
Z_grid = NaN*ones(length(X(:,1)),length(X(1,:)));

v = VideoWriter('ZigZagGrid50PointsWithCavity.avi');
open(v);

File_name = ['.\Data Analysis\Outline_LAGOS.mat'];
Outline = MM*importdata(File_name);

for d=1:length(Sorted_plot_array(:,1))
    [tf, index]=ismember(Sorted_plot_array(d,:),Results_array_shift,'rows');

    Index_x = ceil(index/Number_of_stations(2));

    if mod(index,Station_number_y+1)==0
        Index_y = Station_number_y;
    else
        Index_y = mod(index, Station_number_y+1) ;
    end

    Z_grid(Index_x,Index_y) = Z(Index_x,Index_y);

    Z_grid_int = inpaint_nans(Z_grid,3);

    if d==1
       base = Z_grid_int; 
    end

    if mod(d,1) == 0
%         figure()
%         contour(X,Y,Z_spiral_int-base)
        
        figure()
        hold on
        
        Contour_levels = 0:0.001:2;
        h1 = contour(X,Y,Z_grid_int,Contour_levels);
        h2 = scatter(Sorted_plot_array(1:d,1),Sorted_plot_array(1:d,2),25,'k','o','filled');
        h3 = plot(Outline(:,1),Outline(:,2),'r--');
        caxis([0,2])
        c = colorbar;
        c.Label.String = 'Gravity anomaly (mGal)';
        
        str = sprintf('Contour map after %d gravity measurements, zig-zag',d);
        title(str,'interpreter','latex')
        xlabel('x-coordinate [m]','interpreter','latex')
        ylabel('y-coordinate [m]','interpreter','latex')
        legend([h2,h3],'Measurement points','Cavity outline', 'Location', 'northwest')
        hold off
        
        frame = getframe(gcf);
        for t=1:12
            writeVideo(v,frame);
        end
        
        MovieMatrix(:,round(d/1)) = getframe(gcf);
%         figure()
%         contour(X,Y,Z_spiral)
    end
end   

close(v)

