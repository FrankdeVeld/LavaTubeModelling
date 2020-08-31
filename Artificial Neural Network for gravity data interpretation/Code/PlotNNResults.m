% Outputs = NNTest1_outputs';
% True_output = Labels';
% 
% figure()
% hold on
% grid on
% scatter(True_output(:,1)/140, Outputs(:,1)/140, 'filled')
% scatter(True_output(:,1)/140, True_output(:,1)/140, 'filled')
% title('Result accuracy of the ANN: Normalized Depth (MatLab ANN)')
% xlabel('True Values [m]','interpreter','latex')
% ylabel('ANN output values [m]','interpreter','latex')
% 
% figure() 
% hold on
% grid on
% scatter(True_output(:,2), Outputs(:,2),'filled')
% title('Result accuracy of the ANN: Normalized Shape (MatLab ANN)')
% xlabel('True Values [-]','interpreter','latex')
% ylabel('ANN output values [-]','interpreter','latex')

Outputs_2 = NN_test2_outputs';
True_output = Labels';

figure()
hold on
grid on
scatter(True_output(:,1)/200, Outputs_2(:,1)/200, 'filled')
scatter(True_output(:,1)/200, True_output(:,1)/200, 'filled')
title('Result accuracy of the ANN: Normalized Depth (MatLab ANN)')
xlabel('True Values [m]','interpreter','latex')
ylabel('ANN output values [m]','interpreter','latex')

figure() 
hold on
grid on
scatter(True_output(:,2), Outputs_2(:,2),'filled')
title('Result accuracy of the ANN: Normalized Shape (MatLab ANN)')
xlabel('True Values [-]','interpreter','latex')
ylabel('ANN output values [-]','interpreter','latex')