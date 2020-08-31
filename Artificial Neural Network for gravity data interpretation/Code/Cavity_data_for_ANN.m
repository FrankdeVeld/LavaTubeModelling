Results_40 = importdata('D:\Even Dropbox Stage zooi\Stage\Documenten Lava Tubes\IGMAS\Database\Database Result Files\Results\LAGOS_MM_1_MD_3_CD_40_SSF_4_3D.mat');

Max_g = 0;
g_index = 0;
for a=1:length(Results_40(:,1))
    if Results_40(a,3)>Max_g
        Max_g = Results_40(a,3);
        g_index = a;
    end
end
for a=1:length(Results_40(:,1))
    Results_40(a,3) = Results_40(a,3)/Max_g;
end

% Picking random measurement points, assuming a sphere or vertical cylinder
Features = [];
Labels = [];
Max_x = 0;
for a=1:50
    rand_loc = ceil(length(Results_40(:,1))*rand);
    x=sqrt((Results_40(rand_loc,1)-Results_40(g_index,1))^2+(Results_40(rand_loc,2)-Results_40(g_index,2))^2);
    if x>Max_x
        Max_x=x;
    end
    Features(a,1) = x;
    Features(a,2) = Results_40(rand_loc,3); 
end

for a=1:50
    Features(a,1) = Features(a,1)/Max_x;
    Labels(a,1) = 40;
end
dlmwrite('Features_Cave_40.txt',Features)
save ( 'Features_Cave_40.txt','Features', '-ascii')

% Picking random measurement points, assuming a horizontal cylinder
Features_HC = [];
Labels_HC = [];
Max_x = 0;
for a=1:50
    rand_loc = ceil(length(Results_40(:,1))*rand);
    x=abs(sqrt((Results_40(a,1)+Results_40(a,2))^2)*sin(Results_40(a,2)/Results_40(a,1)-atan(1)));
    if x>Max_x
        Max_x=x;
    end
    Features_HC(a,1) = x;
    Features_HC(a,2) = Results_40(rand_loc,3); 
end

for a=1:50
    Features_HC(a,1) = Features_HC(a,1)/Max_x;
    Labels_HC(a,1) = 40;
end

dlmwrite('Features_HC_Cave_40.txt',Features_HC)
save ( 'Features_HC_Cave_40.txt','Features_HC', '-ascii')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results_80 = importdata('D:\Even Dropbox Stage zooi\Stage\Documenten Lava Tubes\IGMAS\Database\Database Result Files\Results\LAGOS_MM_1_MD_3_CD_40_SSF_4_3D.mat');

Max_g = 0;
g_index = 0;
for a=1:length(Results_80(:,1))
    if Results_80(a,3)>Max_g
        Max_g = Results_80(a,3);
        g_index = a;
    end
end
for a=1:length(Results_80(:,1))
    Results_80(a,3) = Results_80(a,3)/Max_g;
end

% Picking random measurement points, assuming a sphere or vertical cylinder
Features = [];
Labels = [];
for a=1:50
    rand_loc = ceil(length(Results_80(:,1))*rand);
    x=sqrt((Results_80(rand_loc,1)-Results_80(g_index,1))^2+(Results_80(rand_loc,2)-Results_80(g_index,2))^2);
    Features(a,1) = x;
    Features(a,2) = Results_80(rand_loc,3); 
end

for a=1:50
    Features(a,1) = Features(a,1)/Max_x;
    Labels(a,1) = 80;
end
dlmwrite('Features_Cave_80.txt',Features)
save ( 'Features_Cave_80.txt','Features', '-ascii')

% Picking random measurement points, assuming a horizontal cylinder
Features_HC = [];
for a=1:50
    rand_loc = ceil(length(Results_80(:,1))*rand);
    x=abs(sqrt((Results_80(a,1)+Results_80(a,2))^2)*sin(Results_80(a,2)/Results_80(a,1)-atan(1)));
    Features_HC(a,1) = x;
    Features_HC(a,2) = Results_40(rand_loc,3); 
end

for a=1:50
    Features_HC(a,1) = Features_HC(a,1)/Max_x;
end

dlmwrite('Features_HC_Cave_80.txt',Features_HC)
save ( 'Features_HC_Cave_80.txt','Features_HC', '-ascii')