%% Sphere
True_output_Sphere_5pn = Labels_Sphere_5pnoise';

NN_Outputs_SPSP_5pn = NN_Sphere_1_Sphere_5pn_outputs';

subplot(3,3,1)
hold on
grid on
scatter(True_output_Sphere_5pn(:,1), NN_Outputs_SPSP_5pn(:,1), 'filled')
scatter(True_output_Sphere_5pn(:,1), True_output_Sphere_5pn(:,1), 'filled')
title('SP data on SP network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_SPHC_5pn = NN_Sphere_1_HCyl_outputs';

subplot(3,3,2)
hold on
grid on
scatter(True_output_Sphere_5pn(:,1), NN_Outputs_SPHC_5pn(:,1), 'filled')
scatter(True_output_Sphere_5pn(:,1), True_output_Sphere_5pn(:,1), 'filled')
title('HC data on SP network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_SPVC_5pn = NN_Sphere_1_VCyl_outputs';

subplot(3,3,3)
hold on
grid on
scatter(True_output_Sphere_5pn(:,1), NN_Outputs_SPVC_5pn(:,1), 'filled')
scatter(True_output_Sphere_5pn(:,1), True_output_Sphere_5pn(:,1), 'filled')
title('VC data on SP network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

%% HCyl

True_output_HCyl_5pn = Labels_HCyl_5pnoise';

NN_Outputs_HCSP_5pn = NN_HCyl_1_Sphere_5pn_outputs';

subplot(3,3,4)
hold on
grid on
scatter(True_output_HCyl_5pn(:,1), NN_Outputs_HCSP_5pn(:,1), 'filled')
scatter(True_output_HCyl_5pn(:,1), True_output_HCyl_5pn(:,1), 'filled')
title('SP data on HC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_HCHC_5pn = NN_HCyl_1_HCyl_5pn_outputs';

subplot(3,3,5)
hold on
grid on
scatter(True_output_HCyl_5pn(:,1), NN_Outputs_HCHC_5pn(:,1), 'filled')
scatter(True_output_HCyl_5pn(:,1), True_output_HCyl_5pn(:,1), 'filled')
title('HC data on HC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_HCVC_5pn = NN_HCyl_1_VCyl_5pn_outputs';

subplot(3,3,6)
hold on
grid on
scatter(True_output_HCyl_5pn(:,1), NN_Outputs_HCVC_5pn(:,1), 'filled')
scatter(True_output_HCyl_5pn(:,1), True_output_HCyl_5pn(:,1), 'filled')
title('VC data on HC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

%% VCyl

True_output_VCyl_5pn = Labels_VCyl_5pnoise';

NN_Outputs_VCSP_5pn = NN_VCyl_1_Sphere_5pn_outputs';

subplot(3,3,7)
hold on
grid on
scatter(True_output_VCyl_5pn(:,1), NN_Outputs_VCSP_5pn(:,1), 'filled')
scatter(True_output_VCyl_5pn(:,1), True_output_VCyl_5pn(:,1), 'filled')
title('SP data on VC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_VCHC_5pn = NN_VCyl_1_HCyl_5pn_outputs';

subplot(3,3,8)
hold on
grid on
scatter(True_output_VCyl_5pn(:,1), NN_Outputs_VCHC_5pn(:,1), 'filled')
scatter(True_output_VCyl_5pn(:,1), True_output_VCyl_5pn(:,1), 'filled')
title('HC data on VC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

NN_Outputs_VCVC_5pn = NN_VCyl_1_VCyl_5pn_outputs';

subplot(3,3,9)
hold on
grid on
scatter(True_output_VCyl_5pn(:,1), NN_Outputs_VCVC_5pn(:,1), 'filled')
scatter(True_output_VCyl_5pn(:,1), True_output_VCyl_5pn(:,1), 'filled')
title('VC data on VC network')
xlabel('True depth [m]','interpreter','latex')
ylabel('ANN output [m]','interpreter','latex')

