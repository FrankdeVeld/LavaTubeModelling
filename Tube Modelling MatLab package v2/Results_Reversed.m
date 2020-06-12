%% Goal of this script: visualise the obtained results in a clear way
% Authors: Frank de Veld
% Date: June 10th, 2020
% Output: 3D colour graph, contourplot

% This file is much work-in-progress and is really just to quickly look at
% results. Feel fre to change whatever you need to show what you want to be
% shown.

% For 'xml2struct', see https://nl.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
With_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package v2\Verification Tests\Vis_Test_With_Cavity_Results.stations');

With_cavity_results_array = [];
for a=1:length(With_cavity_results.geodata.vertex(1,:))

    Temp_x = With_cavity_results.geodata.vertex{1,a}.Attributes.x;
    Temp_y = With_cavity_results.geodata.vertex{1,a}.Attributes.y;
    
    % If you have calculated sevaral anomalies instead, change 'property'
    % to 'property{1,1}', and the second 1 can be changed to 2, 3, 4, etc.,
    % what you want to see. Most likely, this sequency is the sequency you
    % chose to calculate gz, gz, hgz, etc.
    Temp_result = With_cavity_results.geodata.vertex{1,a}.property.Attributes.value;

    With_cavity_results_array(a,1) = str2num(Temp_x);
    With_cavity_results_array(a,2) = str2num(Temp_y);
    With_cavity_results_array(a,3) = str2num(Temp_result);
end
clear a

No_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package v2\Verification Tests\Vis_Test_No_Cavity_Results.stations');
No_cavity_results_array = [];
for a=1:length(With_cavity_results.geodata.vertex(1,:))

    Temp_x = No_cavity_results.geodata.vertex{1,a}.Attributes.x;
    Temp_y = No_cavity_results.geodata.vertex{1,a}.Attributes.y;
    
    % If you have calculated sevaral anomalies instead, change 'property'
    % to 'property{1,1}', and the second 1 can be changed to 2, 3, 4, etc.,
    % what you want to see. Most likely, this sequency is the sequency you
    % chose to calculate gz, gz, hgz, etc.
    Temp_result = No_cavity_results.geodata.vertex{1,a}.property.Attributes.value;

    No_cavity_results_array(a,1) = str2num(Temp_x);
    No_cavity_results_array(a,2) = str2num(Temp_y);
    No_cavity_results_array(a,3) = str2num(Temp_result);
end
clear a

Reversed_cavity_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package v2\Verification Tests\Vis_Test_Reversed_With_Cavity_Results.stations');

Reversed_cavity_results_array = [];
for a=1:length(Reversed_cavity_results.geodata.vertex(1,:))

    Temp_x = Reversed_cavity_results.geodata.vertex{1,a}.Attributes.x;
    Temp_y = Reversed_cavity_results.geodata.vertex{1,a}.Attributes.y;
    
    % If you have calculated sevaral anomalies instead, change 'property'
    % to 'property{1,1}', and the second 1 can be changed to 2, 3, 4, etc.,
    % what you want to see. Most likely, this sequency is the sequency you
    % chose to calculate gz, gz, hgz, etc.
    Temp_result = Reversed_cavity_results.geodata.vertex{1,a}.property.Attributes.value;

    Reversed_cavity_results_array(a,1) = str2num(Temp_x);
    Reversed_cavity_results_array(a,2) = str2num(Temp_y);
    Reversed_cavity_results_array(a,3) = str2num(Temp_result);
end
clear a

% With the results of both files, we can make a graph of the difference
% between the models with and without cavity. Note that the range is larger
% than the one of the cavity for reducing edge effects. However, since the
% edge effect is present in both models, taking the difference (almost) eliminates
% that.
tri = delaunay(With_cavity_results_array(:,1), With_cavity_results_array(:,2));
figure()
axis tight
trisurf(tri, With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3)) %No_cavity_results_array(:,3)-With_cavity_results_array(:,3)
title('Signal cavity','interpreter','latex')
xlabel('x-coordinate [m]','interpreter','latex')
ylabel('y-coordinate [m]','interpreter','latex')
zlabel('Gravity Difference [mGal]','interpreter','latex')
c = colorbar;
c.Label.String = 'Gravity Difference [mGal]';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 11;

labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [mGal]','Signal cavity' };
ContourPlot(2,With_cavity_results_array(:,1),With_cavity_results_array(:,2),No_cavity_results_array(:,3)-With_cavity_results_array(:,3),labels) %No_cavity_results_array(:,3)-With_cavity_results_array(:,3)

%

% With the results of both files, we can make a graph of the difference
% between the models with and without cavity. Note that the range is larger
% than the one of the cavity for reducing edge effects. However, since the
% edge effect is present in both models, taking the difference (almost) eliminates
% that.
tri = delaunay(With_cavity_results_array(:,1), With_cavity_results_array(:,2));
figure()
axis tight
trisurf(tri, With_cavity_results_array(:,1),With_cavity_results_array(:,2),Reversed_cavity_results_array(:,3)) %No_cavity_results_array(:,3)-With_cavity_results_array(:,3)
title('Signal reversed cavity','interpreter','latex')
xlabel('x-coordinate [m]','interpreter','latex')
ylabel('y-coordinate [m]','interpreter','latex')
zlabel('Gravity Difference [mGal]','interpreter','latex')
c = colorbar;
c.Label.String = 'Gravity Difference [mGal]';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 11;

labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [mGal]','Signal reversed cavity' };
ContourPlot(2,With_cavity_results_array(:,1),With_cavity_results_array(:,2),Reversed_cavity_results_array(:,3),labels) %No_cavity_results_array(:,3)-With_cavity_results_array(:,3)



%% Plotting Functions
function ContourPlot(fignumber,x,y,z,labels)
figure()
xlin = linspace(min(x),max(x),2500);                         
ylin = linspace(min(y),max(y),2500);                         
[X,Y] = meshgrid(xlin,ylin);                                      
Z = griddata(x,y,z,X,Y,'linear');
contourf(X,Y,Z, 'LineColor', 'none' )
hold on
xlabel(labels(1),'interpreter','latex','fontsize',12)
ylabel(labels(2),'interpreter','latex','fontsize',12)
cb=colorbar;
cb.Label.String = labels(3);
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 12;
caxis([min(z) max(z)])
title(labels(4),'interpreter','latex','fontsize',12)
set(gcf,'position',[100,100,800,400])
end

function ContourSubPlot(fignumber,x,y,z,labels)
figure()
xlin = linspace(min(x),max(x),2500);                         
ylin = linspace(min(y),max(y),2500);                         
[X,Y] = meshgrid(xlin,ylin);                                      
Z = griddata(x,y,z,X,Y,'linear');
contourf(X,Y,Z, 'LineColor', 'none' )
hold on
xlabel(labels(1),'interpreter','latex','fontsize',12)
ylabel(labels(2),'interpreter','latex','fontsize',12)
cb=colorbar;
cb.Label.String = labels(3);
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 12;
caxis([min(z) max(z)])
set(gcf,'position',[100,100,600,300])
end
