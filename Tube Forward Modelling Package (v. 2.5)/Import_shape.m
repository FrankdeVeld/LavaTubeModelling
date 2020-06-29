function Cavity_coordinates = Import_shape(Shape, Length, Width, Height, Radius)
    if strcmp(Shape,'Cylinder_H')==true
        [x,z,y] = cylinder(Radius);
        y = y*Length;
        Cavity_coordinates = [x(1,:)',y(1,:)',z(1,:)';x(2,:)',y(2,:)',z(2,:)'];
    elseif strcmp(Shape,'Cylinder_V')==true
        [x,y,z] = cylinder(Radius,20);
        z = z*Height;
        Cavity_coordinates = [x(1,:)',y(1,:)',z(1,:)';x(2,:)',y(2,:)',z(2,:)'];
    elseif strcmp(Shape,'Sphere')==true
        [x,y,z] = sphere;
        x = x * Radius;
        y = y * Radius;
        z = z * Radius;
        
        Cavity_coordinates = [];
        for a = 1:length(x(:,1))
            Temp_coords = [x(a,:)',y(a,:)',z(a,:)'];
            Cavity_coordinates = [Cavity_coordinates;Temp_coords];
        end
    elseif strcmp(Shape,'Beam')==true
        for a=1:8
            if (a == 1 || a == 5)
                Cavity_coordinates(a,1) = 0;
                Cavity_coordinates(a,2) = 0;    
            elseif (a==2 || a == 6)
                Cavity_coordinates(a,1) = Length;
                Cavity_coordinates(a,2) = 0;
            elseif (a==3 || a == 7)
                Cavity_coordinates(a,1) = 0;
                Cavity_coordinates(a,2) = Width;
            elseif (a==4 || a == 8)
                Cavity_coordinates(a,1) = Length;
                Cavity_coordinates(a,2) = Width;
            end
            if a<5
                Cavity_coordinates(a,3) = 0;
            else
                Cavity_coordinates(a,3) = Height;
            end
        end
    end
end