With_cavity_results_sh3 = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package\Sizing Tests and Results\sh3_With_Cavity.stations');

With_cavity_results_array_sh3 = [];
for a=1:length(With_cavity_results_sh3.geodata.vertex(1,:))

    Temp_x = With_cavity_results_sh3.geodata.vertex{1,a}.Attributes.x;
    Temp_y = With_cavity_results_sh3.geodata.vertex{1,a}.Attributes.y;

    Temp_result = With_cavity_results_sh3.geodata.vertex{1,a}.property.Attributes.value;

    With_cavity_results_array_sh3(a,1) = str2num(Temp_x);
    With_cavity_results_array_sh3(a,2) = str2num(Temp_y);
    With_cavity_results_array_sh3(a,3) = str2num(Temp_result);
end
clear a

No_cavity_results_sh3 = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package\Sizing Tests and Results\sh3_No_Cavity.stations');

No_cavity_results_array_sh3 = [];
for a=1:length(No_cavity_results_sh3.geodata.vertex(1,:))

    Temp_x = No_cavity_results_sh3.geodata.vertex{1,a}.Attributes.x;
    Temp_y = No_cavity_results_sh3.geodata.vertex{1,a}.Attributes.y;

    Temp_result = No_cavity_results_sh3.geodata.vertex{1,a}.property.Attributes.value;

    No_cavity_results_array_sh3(a,1) = str2num(Temp_x);
    No_cavity_results_array_sh3(a,2) = str2num(Temp_y);
    No_cavity_results_array_sh3(a,3) = str2num(Temp_result);
end
clear a

tri = delaunay(With_cavity_results_array_sh3(:,1), With_cavity_results_array_sh3(:,2));
figure()
axis tight
trisurf(tri, With_cavity_results_array_sh3(:,1),With_cavity_results_array_sh3(:,2),No_cavity_results_array_sh3(:,3)-With_cavity_results_array_sh3(:,3))
title('Difference between the signals of models with and without a cavity, Stretching Horizontal 3','interpreter','latex')
xlabel('x-coordinate [m]','interpreter','latex')
ylabel('y-coordinate [m]','interpreter','latex')
zlabel('Gravity Difference [mGal]','interpreter','latex')
c = colorbar;
c.Label.String = 'Gravity Difference [mGal]';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 11;

labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [mGal]','Difference between the signals of models with and without a cavity, Stretching Horizontal 3' };
ContourPlot(2,With_cavity_results_array_sh3(:,1),With_cavity_results_array_sh3(:,2),No_cavity_results_array_sh3(:,3)-With_cavity_results_array_sh3(:,3),labels)



With_cavity_results_Def = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package\Sizing Tests and Results\Def_With_Cavity.stations');

With_cavity_results_array_Def = [];
for a=1:length(With_cavity_results_Def.geodata.vertex(1,:))

    Temp_x = With_cavity_results_Def.geodata.vertex{1,a}.Attributes.x;
    Temp_y = With_cavity_results_Def.geodata.vertex{1,a}.Attributes.y;

    Temp_result = With_cavity_results_Def.geodata.vertex{1,a}.property.Attributes.value;

    With_cavity_results_array_Def(a,1) = str2num(Temp_x);
    With_cavity_results_array_Def(a,2) = str2num(Temp_y);
    With_cavity_results_array_Def(a,3) = str2num(Temp_result);
end
clear a

No_cavity_results_Def = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package\Sizing Tests and Results\Def_No_Cavity.stations');

No_cavity_results_array_Def = [];
for a=1:length(No_cavity_results_Def.geodata.vertex(1,:))

    Temp_x = No_cavity_results_Def.geodata.vertex{1,a}.Attributes.x;
    Temp_y = No_cavity_results_Def.geodata.vertex{1,a}.Attributes.y;

    Temp_result = No_cavity_results_Def.geodata.vertex{1,a}.property.Attributes.value;

    No_cavity_results_array_Def(a,1) = str2num(Temp_x);
    No_cavity_results_array_Def(a,2) = str2num(Temp_y);
    No_cavity_results_array_Def(a,3) = str2num(Temp_result);
end
clear a


tri = delaunay(With_cavity_results_array_Def(:,1), With_cavity_results_array_Def(:,2));
figure()
axis tight
trisurf(tri, With_cavity_results_array_Def(:,1),With_cavity_results_array_Def(:,2),No_cavity_results_array_Def(:,3)-With_cavity_results_array_Def(:,3))
title('Difference between the signals of models with and without a cavity, default','interpreter','latex')
xlabel('x-coordinate [m]','interpreter','latex')
ylabel('y-coordinate [m]','interpreter','latex')
zlabel('Gravity Difference [mGal]','interpreter','latex')
c = colorbar;
c.Label.String = 'Gravity Difference [mGal]';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 11;

labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [mGal]','Difference between the signals of models with and without a cavity, default' };
ContourPlot(4,With_cavity_results_array_Def(:,1),With_cavity_results_array_Def(:,2),No_cavity_results_array_Def(:,3)-With_cavity_results_array_Def(:,3),labels)

function ContourPlot(fignumber,x,y,z,labels)
figure(fignumber)
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
figure(fignumber)
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