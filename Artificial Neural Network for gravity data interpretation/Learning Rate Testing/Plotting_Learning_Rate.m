% Plotting for the learning rate investigation
% Constant: Data points 100x15, epochs 10, hidden neurons 64
% Variable: initial learning rate

Learning_rates = [0.1,0.25,0.5,1];
Av_tot_spread_test = [0.20064648287943543, 0.2637419230171778, 0.42387781882684983, 0.751093707294318];
Av_std_tot_spread_test = [0.03117120847025674, 0.036315773413089086, 0.052974678872124, 0.08742165746710026];
Av_tot_spread_train = [0.20865563877290122, 0.2740693696122882, 0.43397732755558305, 0.7448629124014334];
Av_std_tot_spread_train = [0.03288671394725525, 0.0376847964907328, 0.05060864139842198, 0.08128232675861599];
Final_loss = [0.038721739537134196, 0.036751626518042875, 0.04269230204119229, 0.053612449145816345];
Final_loss_val = [0.04168517032317177, 0.037002399438320545, 0.04756331021598333, 0.07057534183649174];

figure()
hold on
grid on
scatter(Learning_rates,Av_tot_spread_test,30,'r','d','filled')
scatter(Learning_rates,Av_tot_spread_train,30,'b','d','filled')
scatter(Learning_rates,Av_std_tot_spread_test,30,'r','filled')
scatter(Learning_rates,Av_std_tot_spread_train,30,'b','filled')
plot(Learning_rates,(Av_tot_spread_test+Av_tot_spread_train)./2,'k--')
plot(Learning_rates,(Av_std_tot_spread_test+Av_std_tot_spread_train)./2,'k--')
title('Performance metrics, spreading, for learning rate variation')
xlabel('Learning rate (-)','interpreter','latex')
ylabel('Average spread / Average STD spread (-)','interpreter','latex')
legend('Average spread, testing','Average spread, training','Average STD spread, testing','Average STD spread, training')

figure()
hold on
grid on
scatter(Learning_rates,Final_loss,30,'r','d','filled')
scatter(Learning_rates,Final_loss_val,30,'b','d','filled')
plot(Learning_rates,(Final_loss+Final_loss_val)./2,'k--')
title('Performance metrics, loss, for learning rate variation')
xlabel('Learning rate (-)','interpreter','latex')
ylabel('Final loss value (-)','interpreter','latex')
legend('Final loss value - training','Final loss value - validating','Average trend')

