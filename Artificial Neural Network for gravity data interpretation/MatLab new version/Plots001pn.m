%% Sphere
True_output_Sphere_001pn = Labels_Sphere_001pnoise';

NN_Outputs_SPSP_001pn = NN_Sphere_1_Sphere_001pn_outputs';

subplot(3,3,1)
hold on
grid on
scatter(True_output_Sphere_001pn(:,1), NN_Outputs_SPSP_001pn(:,1), 'filled')
scatter(True_output_Sphere_001pn(:,1), True_output_Sphere_001pn(:,1), 'filled')
title('SP data on SP network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_SPHC_001pn = NN_Sphere_1_HCyl_outputs';

subplot(3,3,2)
hold on
grid on
scatter(True_output_Sphere_001pn(:,1), NN_Outputs_SPHC_001pn(:,1), 'filled')
scatter(True_output_Sphere_001pn(:,1), True_output_Sphere_001pn(:,1), 'filled')
title('HC data on SP network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_SPVC_001pn = NN_Sphere_1_VCyl_outputs';

subplot(3,3,3)
hold on
grid on
scatter(True_output_Sphere_001pn(:,1), NN_Outputs_SPVC_001pn(:,1), 'filled')
scatter(True_output_Sphere_001pn(:,1), True_output_Sphere_001pn(:,1), 'filled')
title('VC data on SP network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

%% HCyl

True_output_HCyl_001pn = Labels_HCyl_001pnoise';

NN_Outputs_HCSP_001pn = NN_HCyl_1_Sphere_001pn_outputs';

subplot(3,3,4)
hold on
grid on
scatter(True_output_HCyl_001pn(:,1), NN_Outputs_HCSP_001pn(:,1), 'filled')
scatter(True_output_HCyl_001pn(:,1), True_output_HCyl_001pn(:,1), 'filled')
title('SP data on HC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_HCHC_001pn = NN_HCyl_1_HCyl_001pn_outputs';

subplot(3,3,5)
hold on
grid on
scatter(True_output_HCyl_001pn(:,1), NN_Outputs_HCHC_001pn(:,1), 'filled')
scatter(True_output_HCyl_001pn(:,1), True_output_HCyl_001pn(:,1), 'filled')
title('HC data on HC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_HCVC_001pn = NN_HCyl_1_VCyl_001pn_outputs';

subplot(3,3,6)
hold on
grid on
scatter(True_output_HCyl_001pn(:,1), NN_Outputs_HCVC_001pn(:,1), 'filled')
scatter(True_output_HCyl_001pn(:,1), True_output_HCyl_001pn(:,1), 'filled')
title('VC data on HC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

%% VCyl

True_output_VCyl_001pn = Labels_VCyl_001pnoise';

NN_Outputs_VCSP_001pn = NN_VCyl_1_Sphere_001pn_outputs';

subplot(3,3,7)
hold on
grid on
scatter(True_output_VCyl_001pn(:,1), NN_Outputs_VCSP_001pn(:,1), 'filled')
scatter(True_output_VCyl_001pn(:,1), True_output_VCyl_001pn(:,1), 'filled')
title('SP data on VC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_VCHC_001pn = NN_VCyl_1_HCyl_001pn_outputs';

subplot(3,3,8)
hold on
grid on
scatter(True_output_VCyl_001pn(:,1), NN_Outputs_VCHC_001pn(:,1), 'filled')
scatter(True_output_VCyl_001pn(:,1), True_output_VCyl_001pn(:,1), 'filled')
title('HC data on VC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_VCVC_001pn = NN_VCyl_1_VCyl_001pn_outputs';

subplot(3,3,9)
hold on
grid on
scatter(True_output_VCyl_001pn(:,1), NN_Outputs_VCVC_001pn(:,1), 'filled')
scatter(True_output_VCyl_001pn(:,1), True_output_VCyl_001pn(:,1), 'filled')
title('VC data on VC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

