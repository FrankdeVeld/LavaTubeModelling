% Shape = 'Sphere';
% Shape = 'Cylinder_H';
% Shape = 'Beam';
Shape = 'Cylinder_V';

if strcmp('Sphere',Shape)==true
    Bare_IGMAS_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package v2\Verification Tests\Bare_Sphere_verification_test_results.stations');
    
    Theoretical_bare_results_array = [];
    for a=1:length(Bare_IGMAS_results.geodata.vertex(1,:))

        x = str2num(Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.x);
        y = str2num(Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.y);
        
        r = sqrt(x^2+y^2);
        rho = 2.5*10^3;
        R = 10;
        z = 10+R;
        G = 6.67408e-11;
        
        % http://www.cas.usf.edu/~cconnor/pot_fields_lectures/Lecture4_gravity.pdf
        gz = 4*pi*G*rho*R^3/3*(z)/((r^2+z^2)^(3/2))*10^5;

        Theoretical_bare_results_array(a,1) = x;
        Theoretical_bare_results_array(a,2) = y;
        Theoretical_bare_results_array(a,3) = gz;
    end
    
    Bare_IGMAS_results_array = [];
    for a=1:length(Bare_IGMAS_results.geodata.vertex(1,:))

        Temp_x = Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.y;
        
        Temp_result = Bare_IGMAS_results.geodata.vertex{1,a}.property.Attributes.value;

        Bare_IGMAS_results_array(a,1) = str2num(Temp_x);
        Bare_IGMAS_results_array(a,2) = str2num(Temp_y);
        Bare_IGMAS_results_array(a,3) = str2num(Temp_result);
    end

elseif strcmp('Cylinder_H',Shape)==true
    Bare_IGMAS_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package v2\Verification Tests\Bare_Cylinder_H_verification_test_results.stations');
    
    Theoretical_bare_results_array = [];
    for a=1:length(Bare_IGMAS_results.geodata.vertex(1,:))

        x = str2num(Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.x);
        y = str2num(Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.y);
        
        rho = 2.5*10^3;
        R = 15;
        z = 10+R;
        L = 100;
        G = 6.67408e-11;
        
        % https://www.wolframalpha.com/input/?i=Integrate%5B1%2F%28z%5E2%2B%28l-x%29%5E2%2By%5E2%29%5E%283%2F2%29%2Cl%5D
        gz = G*pi*R^2*rho*z*((L-y)/((z^2+x^2)*sqrt((L-y)^2+x^2+z^2)) + y/((z^2+x^2)*sqrt(z^2 +y^2+x^2))) *10^5;

        Theoretical_bare_results_array(a,1) = x;
        Theoretical_bare_results_array(a,2) = y;
        Theoretical_bare_results_array(a,3) = gz;
    end
    
    Bare_IGMAS_results_array = [];
    for a=1:length(Bare_IGMAS_results.geodata.vertex(1,:))

        Temp_x = Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.y;
        
        Temp_result = Bare_IGMAS_results.geodata.vertex{1,a}.property.Attributes.value;

        Bare_IGMAS_results_array(a,1) = str2num(Temp_x);
        Bare_IGMAS_results_array(a,2) = str2num(Temp_y);
        Bare_IGMAS_results_array(a,3) = str2num(Temp_result);
    end
elseif strcmp('Cylinder_V',Shape)==true
    Bare_IGMAS_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package v2\Verification Tests\Bare_Cylinder_V_verification_test_results.stations');
    
    Theoretical_bare_results_array = [];
    for a=1:length(Bare_IGMAS_results.geodata.vertex(1,:))

        x = str2num(Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.x);
        y = str2num(Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.y);
        
        r = sqrt(x^2+y^2);
        rho = 2.5*10^3;
        R = 15;
        z = 10;
        L = 100;
        G = 6.67408e-11;
        % http://www.cas.usf.edu/~cconnor/pot_fields_lectures/Lecture4_gravity.pdf
        gz = G*pi*R^2*rho*(1/sqrt(r^2+z^2) - 1/sqrt(r^2 + (z+L)^2)) *10^5;

        Theoretical_bare_results_array(a,1) = x;
        Theoretical_bare_results_array(a,2) = y;
        Theoretical_bare_results_array(a,3) = gz;
    end
    
    
    Bare_IGMAS_results_array = [];
    for a=1:length(Bare_IGMAS_results.geodata.vertex(1,:))

        Temp_x = Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.y;
        
        Temp_result = Bare_IGMAS_results.geodata.vertex{1,a}.property.Attributes.value;

        Bare_IGMAS_results_array(a,1) = str2num(Temp_x);
        Bare_IGMAS_results_array(a,2) = str2num(Temp_y);
        Bare_IGMAS_results_array(a,3) = str2num(Temp_result);
    end
    
    
elseif strcmp('Beam',Shape)==true
    
    Bare_IGMAS_results = xml2struct('C:\Users\frank\Dropbox\Studie\Stage\Documenten Lava Tubes\IGMAS\Tube Modelling MatLab package v2\Verification Tests\Bare_Beam_verification_test_results.stations');
    
    Theoretical_bare_results_array = [];
    for a=1:length(Bare_IGMAS_results.geodata.vertex(1,:))

        x = str2num(Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.x);
        y = str2num(Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.y);
        
        r = sqrt(x^2+y^2);
        rho = 2.5*10^3;
        R = 10;
        z = 10+R;
        G = 6.67408e-11;
        
        % This one can't be analytically solved, or not quickly in any case
        gz = 4*pi*G*rho*R^3/3*(z)/((r^2+z^2)^(3/2))*10^5;

        Theoretical_bare_results_array(a,1) = x;
        Theoretical_bare_results_array(a,2) = y;
        Theoretical_bare_results_array(a,3) = gz;
    end
    
    
    Bare_IGMAS_results_array = [];
    for a=1:length(Bare_IGMAS_results.geodata.vertex(1,:))

        Temp_x = Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.x;
        Temp_y = Bare_IGMAS_results.geodata.vertex{1,a}.Attributes.y;
        
        Temp_result = Bare_IGMAS_results.geodata.vertex{1,a}.property.Attributes.value;

        Bare_IGMAS_results_array(a,1) = str2num(Temp_x);
        Bare_IGMAS_results_array(a,2) = str2num(Temp_y);
        Bare_IGMAS_results_array(a,3) = str2num(Temp_result);
    end
end

%% IGMAS results figures

tri = delaunay(Bare_IGMAS_results_array(:,1), Bare_IGMAS_results_array(:,2));
figure()
axis tight
trisurf(tri, Bare_IGMAS_results_array(:,1),Bare_IGMAS_results_array(:,2),Bare_IGMAS_results_array(:,3))
title('IGMAS gravity signal','interpreter','latex')
xlabel('x-coordinate [m]','interpreter','latex')
ylabel('y-coordinate [m]','interpreter','latex')
zlabel('Gravity Difference [mGal]','interpreter','latex')
c = colorbar;
c.Label.String = 'Gravity Difference [mGal]';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 11;

labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [mGal]','IGMAS gravity signal' };
ContourPlot(2,Bare_IGMAS_results_array(:,1),Bare_IGMAS_results_array(:,2),Bare_IGMAS_results_array(:,3),labels) 

%% Theoretical results figures
tri = delaunay( Theoretical_bare_results_array(:,1),  Theoretical_bare_results_array(:,2));
figure()
axis tight
trisurf(tri,  Theoretical_bare_results_array(:,1), Theoretical_bare_results_array(:,2), Theoretical_bare_results_array(:,3))
title('Theoretical gravity signal','interpreter','latex')
xlabel('x-coordinate [m]','interpreter','latex')
ylabel('y-coordinate [m]','interpreter','latex')
zlabel('Gravity Difference [mGal]','interpreter','latex')
c = colorbar;
c.Label.String = 'Gravity Difference [mGal]';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 11;

labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [mGal]','Theoretical gravity signal' };
ContourPlot(2, Theoretical_bare_results_array(:,1), Theoretical_bare_results_array(:,2), Theoretical_bare_results_array(:,3),labels) 

%% Comparison figures
tri = delaunay( Theoretical_bare_results_array(:,1),  Theoretical_bare_results_array(:,2));
figure()
axis tight
trisurf(tri,  Theoretical_bare_results_array(:,1), Theoretical_bare_results_array(:,2),  (Theoretical_bare_results_array(:,3)-Bare_IGMAS_results_array(:,3))./Theoretical_bare_results_array(:,3)*100)
title('Difference between signals','interpreter','latex')
xlabel('x-coordinate [m]','interpreter','latex')
ylabel('y-coordinate [m]','interpreter','latex')
zlabel('Gravity Difference [% of theoretical]','interpreter','latex')
c = colorbar;
c.Label.String = 'Gravity Difference [mGal]';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 11;

labels = {'x-coordinate [m]','y-coordinate [m]','Gravity Difference [% of theoretical]','Difference between signals' };
ContourPlot(2, Theoretical_bare_results_array(:,1), Theoretical_bare_results_array(:,2), (Theoretical_bare_results_array(:,3)-Bare_IGMAS_results_array(:,3))./Theoretical_bare_results_array(:,3)*100,labels) 


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