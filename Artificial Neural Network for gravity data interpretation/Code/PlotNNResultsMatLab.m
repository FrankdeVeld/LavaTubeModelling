Outputs = NNTest1_outputs';
True_output = Labels';

figure()
grid on
scatter(True_output(:,1)/140, Outputs(:,1)/140)
title('Performance metrics, spreading, for data set size variation')
xlabel('Data set size (-)','interpreter','latex')
ylabel('Average spread / Average STD spread (-)','interpreter','latex')

figure() 
grid on
scatter(True_output(:,2), Outputs(:,2))