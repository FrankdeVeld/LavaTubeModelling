% Plotting for the learning rate investigation
% Constant: Data points 100x15, epochs 10, hidden neurons 64
% Variable: initial learning rate

Neuron_number = [2,4,8,16,32,64,128,256,512];
Av_tot_spread = [0.7921511325862359, 0.29715123983069913, 0.2538606081782028, 0.2289873344658262, 0.249174862834, 0.2756846507042808, 0.3475819546558902, 0.36988007280508417, 0.6258931362980307];
Av_std_tot_spread = [0.00047213569005108587, 0.03818755134404724, 0.034037214719353656, 0.029041124644745085, 0.03387426252, 0.037004525914001316, 0.046652157223541324, 0.05166186076158347, 0.0811288646507505];
Final_loss = [0.1719474537332125, 0.04088455881416497, 0.0357879186766016, 0.03614229, 0.03516868, 0.03388486, 0.0347998, 0.03482, 0.0690284493246276];
Final_loss_val = [0.2438388167772029, 0.031160693137729022, 0.04064447143373449, 0.03614229, 0.04364276, 0.03826277, 0.03518061, 0.03516895, 0.08145523443886923];

figure()
hold on
grid on
scatter(Neuron_number,Av_tot_spread,30,'r','d','filled')
scatter(Neuron_number,Av_std_tot_spread,30,'b','filled')
title('Performance metrics, spreading, for learning rate variation')
set(gca,'xscale','log')
xlabel('Number of neurons (-)','interpreter','latex')
ylabel('Average spread / Average STD spread (-)','interpreter','latex')
legend('Average spread','Average STD spread')

figure()
hold on
grid on
scatter(Neuron_number,Final_loss,30,'r','d','filled')
scatter(Neuron_number,Final_loss_val,30,'b','d','filled')
plot(Neuron_number,(Final_loss+Final_loss_val)./2,'k--')
set(gca,'xscale','log')
title('Performance metrics, loss, for learning rate variation')
xlabel('Number of neurons (-)','interpreter','latex')
ylabel('Final loss value (-)','interpreter','latex')
legend('Final loss value - training','Final loss value - validating','Average trend')