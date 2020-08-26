% Plotting for the learning rate investigation
% Constant: Data points 100x15, epochs 10, hidden neurons 64
% Variable: initial learning rate

Epoch_number = [5,10,20,50,100];
Av_tot_spread_test = [0.29801600727718297, 0.2637419230171778, 0.29364985526331594, 0.4836616737951146, 0.821287422063041];
Av_std_tot_spread_test = [0.04598788285742021, 0.036315773413089086, 0.039575348752182424, 0.05914778762560924, 0.10373308718982929];
Av_tot_spread_train = [0.3092055857161695, 0.2740693696122882, 0.3041746104963311, 0.48159281230719525, 0.8093820592476414 ];
Av_std_tot_spread_train = [0.048128360118097166, 0.0376847964907328, 0.03989345053619406, 0.053117126089417356, 0.09328933830482797];
Final_loss = [0.056420072631288946, 0.036751626518042875, 0.035671163046185575, 0.045957951330233345, 0.04384665059130419];
Final_loss_val = [0.05657566614297914, 0.037002399438320545, 0.041396167079469436, 0.06130591463936289, 0.04700898802375925];

figure()
hold on
grid on
scatter(Epoch_number,Av_tot_spread_test,30,'r','d','filled')
scatter(Epoch_number,Av_tot_spread_train,30,'b','d','filled')
scatter(Epoch_number,Av_std_tot_spread_test,30,'r','filled')
scatter(Epoch_number,Av_std_tot_spread_train,30,'b','filled')
plot(Epoch_number,(Av_tot_spread_test+Av_tot_spread_train)./2,'k--')
plot(Epoch_number,(Av_std_tot_spread_test+Av_std_tot_spread_train)./2,'k--')
title('Performance metrics, spreading, for learning rate variation')
xlabel('Number of epochs (-)','interpreter','latex')
ylabel('Average spread / Average STD spread (-)','interpreter','latex')
legend('Average spread, testing','Average spread, training','Average STD spread, testing','Average STD spread, training')

figure()
hold on
grid on
scatter(Epoch_number,Final_loss,30,'r','d','filled')
scatter(Epoch_number,Final_loss_val,30,'b','d','filled')
plot(Epoch_number,(Final_loss+Final_loss_val)./2,'k--')
title('Performance metrics, loss, for learning rate variation')
xlabel('Number of epochs (-)','interpreter','latex')
ylabel('Final loss value (-)','interpreter','latex')
legend('Final loss value - training','Final loss value - validating','Average trend')