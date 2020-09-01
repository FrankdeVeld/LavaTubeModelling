% Author: Frank de Veld
% Code for generating input data for the MatLab ANN
% Version: 1.0
% Date: August 25th, 2020

Features = [];
Labels = [];

for a=1:5000
    z = abs(normrnd(0,50));
    q= ceil(3*rand)/2;
    for b=1:15
        x = abs(normrnd(0,50));
        g = (z^2/(x^2+z^2))^q;
        Features((a-1)*15+b,1) = x;
        Features((a-1)*15+b,2) = g;
        Labels((a-1)*15+b,1) = z;
        Labels((a-1)*15+b,2) = q;
    end
end

Features = Features';
Labels = Labels';