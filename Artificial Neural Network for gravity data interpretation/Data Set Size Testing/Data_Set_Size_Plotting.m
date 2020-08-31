% Plotting for the data set investigation
% Constant: Learning rate 0.25, epochs 10, hidden neurons 64
% Variable: data set size

Data_set_size = [10,50,100,1000,10000];
Av_tot_spread_test = [0.4671557561671985, 0.37301381613566376, 0.2637419230171778, 0.4401521495912329, 2.4334072021193336];
Av_std_tot_spread_test = [0.055098830862402484, 0.06297294146743823, 0.036315773413089086, 0.0497334521005738, 0.30614249656379];
Av_tot_spread_train = [0.4450492709280156, 0.37155301779143446, 0.2740693696122882, 0.4495855079197263, 2.4450386833719127];
Av_std_tot_spread_train = [0.03460992566452257, 0.054639314601349576, 0.0376847964907328, 0.05801886143756689, 0.3093562539364917];
Final_loss = [0.2034260192654926, 0.057399576300880245, 0.036751626518042875, 0.05333092629741505, 0.033581884603366356];
Final_loss_val = [0.1481759466916156, 0.04167237535725126, 0.037002399438320545, 0.06228338957767713,0.033760935303485354];

figure()
hold on
grid on
scatter(Data_set_size,Av_tot_spread_test,30,'r','d','filled')
scatter(Data_set_size,Av_tot_spread_train,30,'b','d','filled')
scatter(Data_set_size,Av_std_tot_spread_test,30,'r','filled')
scatter(Data_set_size,Av_std_tot_spread_train,30,'b','filled')
plot(Data_set_size,(Av_tot_spread_test+Av_tot_spread_train)./2,'k--')
plot(Data_set_size,(Av_std_tot_spread_test+Av_std_tot_spread_train)./2,'k--')
set(gca,'xscale','log')
title('Performance metrics, spreading, for data set size variation')
xlabel('Data set size (-)','interpreter','latex')
ylabel('Average spread / Average STD spread (-)','interpreter','latex')
legend('Average spread, testing','Average spread, training','Average STD spread, testing','Average STD spread, training')

figure()
hold on
grid on
scatter(Data_set_size,Final_loss,30,'r','d','filled')
scatter(Data_set_size,Final_loss_val,30,'b','d','filled')
plot(Data_set_size,(Final_loss+Final_loss_val)./2,'k--')
set(gca,'xscale','log')
title('Performance metrics, loss, for data set size variation')
xlabel('Data set size (-)','interpreter','latex')
ylabel('Final loss value (-)','interpreter','latex')
legend('Final loss value - training','Final loss value - validating','Average trend')