%% Sphere
True_output_Sphere = Labels_Sphere';

NN_Outputs_SPSP = NN_Sphere_1_Sphere_outputs';

subplot(3,3,1)
hold on
grid on
scatter(True_output_Sphere(:,1), NN_Outputs_SPSP(:,1), 'filled')
scatter(True_output_Sphere(:,1), True_output_Sphere(:,1), 'filled')
title('SP data on SP network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_SPHC = NN_Sphere_1_HCyl_outputs';

subplot(3,3,2)
hold on
grid on
scatter(True_output_Sphere(:,1), NN_Outputs_SPHC(:,1), 'filled')
scatter(True_output_Sphere(:,1), True_output_Sphere(:,1), 'filled')
title('HC data on SP network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_SPVC = NN_Sphere_1_VCyl_outputs';

subplot(3,3,3)
hold on
grid on
scatter(True_output_Sphere(:,1), NN_Outputs_SPVC(:,1), 'filled')
scatter(True_output_Sphere(:,1), True_output_Sphere(:,1), 'filled')
title('VC data on SP network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

%% HCyl

True_output_HCyl = Labels_HCyl';

NN_Outputs_HCSP = NN_HCyl_1_Sphere_outputs';

subplot(3,3,4)
hold on
grid on
scatter(True_output_HCyl(:,1), NN_Outputs_HCSP(:,1), 'filled')
scatter(True_output_HCyl(:,1), True_output_HCyl(:,1), 'filled')
title('SP data on HC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_HCHC = NN_HCyl_1_HCyl_outputs';

subplot(3,3,5)
hold on
grid on
scatter(True_output_HCyl(:,1), NN_Outputs_HCHC(:,1), 'filled')
scatter(True_output_HCyl(:,1), True_output_HCyl(:,1), 'filled')
title('HC data on HC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_HCVC = NN_HCyl_1_VCyl_outputs';

subplot(3,3,6)
hold on
grid on
scatter(True_output_HCyl(:,1), NN_Outputs_HCVC(:,1), 'filled')
scatter(True_output_HCyl(:,1), True_output_HCyl(:,1), 'filled')
title('VC data on HC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

%% VCyl

True_output_VCyl = Labels_VCyl';

NN_Outputs_VCSP = NN_VCyl_1_Sphere_outputs';

subplot(3,3,7)
hold on
grid on
scatter(True_output_VCyl(:,1), NN_Outputs_VCSP(:,1), 'filled')
scatter(True_output_VCyl(:,1), True_output_VCyl(:,1), 'filled')
title('SP data on VC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_VCHC = NN_VCyl_1_HCyl_outputs';

subplot(3,3,8)
hold on
grid on
scatter(True_output_VCyl(:,1), NN_Outputs_VCHC(:,1), 'filled')
scatter(True_output_VCyl(:,1), True_output_VCyl(:,1), 'filled')
title('HC data on VC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_VCVC = NN_VCyl_1_VCyl_outputs';

subplot(3,3,9)
hold on
grid on
scatter(True_output_VCyl(:,1), NN_Outputs_VCVC(:,1), 'filled')
scatter(True_output_VCyl(:,1), True_output_VCyl(:,1), 'filled')
title('VC data on VC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')




