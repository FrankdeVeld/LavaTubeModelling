%% Goal of this script: visualise the obtained results in a clear way
% Authors: Frank de Veld
% Date: June 15th, 2020
% Output: 3D colour graph, contourplot

% This file is much work-in-progress and is really just to quickly look at
% results. Feel free to change whatever you need to show what you want to be
% shown.

% For 'xml2struct', see https://nl.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
Rresults = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package v2.5\xxx_results.stations');

Results_array = [];
for a=1:length(Rresults.geodata.vertex(1,:))

    Temp_x = Rresults.geodata.vertex{1,a}.Attributes.x;
    Temp_y = Rresults.geodata.vertex{1,a}.Attributes.y;
    
    % If you have calculated sevaral anomalies instead, change 'property'
    % to 'property{1,1}', and the second 1 can be changed to 2, 3, 4, etc.,
    % what you want to see. Most likely, this sequency is the sequency you
    % chose to calculate gz, gz, hgz, etc.
    Temp_result = Rresults.geodata.vertex{1,a}.property.Attributes.value;

    Results_array(a,1) = str2num(Temp_x);
    Results_array(a,2) = str2num(Temp_y);
    Results_array(a,3) = str2num(Temp_result);
end
clear a

% With the results, we can make a contour graph. Note that the range is larger
% than the one of the cavity for reducing edge effects. 
tri = delaunay(Results_array(:,1), Results_array(:,2));
figure()
axis tight
trisurf(tri, Results_array(:,1),Results_array(:,2),Results_array(:,3))
title('Signal with cavity - signal without cavity','interpreter','latex')
xlabel('x-coordinate [m]','interpreter','latex')
ylabel('y-coordinate [m]','interpreter','latex')
zlabel('Gravity Difference [mGal]','interpreter','latex')
c = colorbar;
c.Label.String = 'Gravity Difference [mGal]';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 11;

labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [mGal]','Signal with cavity - signal without cavity' };
ContourPlot(2,Results_array(:,1),Results_array(:,2),Results_array(:,3),labels) 
