figure
subplot(2,2,1)
semilogx(Tsw,p_hinSigma(1,:),'kx','MarkerSize',10,'LineWidth',2)
hold on
semilogx(Tsw,p_hinSigma(2,:),'rx','MarkerSize',10,'LineWidth',2)

semilogx(Tsw,p_hinSigma(3,:),'ko','MarkerSize',10,'LineWidth',2)
semilogx(Tsw,p_hinSigma(4,:),'ro','MarkerSize',10,'LineWidth',2)

grid on 

% yLim = ylim;
% yLim(1) = 0;
% ylim(yLim)

legend('N Pi-lump','N Pi-lump with unsteady fric.','fMOC','fMOC with unsteady fric.')
ylabel('Pressure, std. dev. (MPa)')
xlabel('Switching period (s)')


subplot(2,2,2)
semilogx(Tsw,p_houtinSigma(1,:),'kx','MarkerSize',10,'LineWidth',2)
hold on
semilogx(Tsw,p_houtinSigma(2,:),'rx','MarkerSize',10,'LineWidth',2)

semilogx(Tsw,p_houtinSigma(3,:),'ko','MarkerSize',10,'LineWidth',2)
semilogx(Tsw,p_houtinSigma(4,:),'ro','MarkerSize',10,'LineWidth',2)

grid on 

% yLim = ylim;
% yLim(1) = 0;
% ylim(yLim)

legend('N Pi-lump','N Pi-lump with unsteady fric.','fMOC','fMOC with unsteady fric.')
ylabel('Pressure, std. dev. (MPa)')
xlabel('Switching period (s)')

subplot(2,2,3)
semilogx(Tsw,p_hinMu(1,:),'kx','MarkerSize',10,'LineWidth',2)
hold on
semilogx(Tsw,p_hinMu(2,:),'rx','MarkerSize',10,'LineWidth',2)

semilogx(Tsw,p_hinMu(3,:),'ko','MarkerSize',10,'LineWidth',2)
semilogx(Tsw,p_hinMu(4,:),'ro','MarkerSize',10,'LineWidth',2)

grid on 

% yLim = ylim;
% yLim(1) = 0;
% ylim(yLim)

legend('N Pi-lump','N Pi-lump with unsteady fric.','fMOC','fMOC with unsteady fric.')
ylabel('Pressure, mean (MPa)')
xlabel('Switching period (s)')


subplot(2,2,4)
semilogx(Tsw,p_houtinMu(1,:),'kx','MarkerSize',10,'LineWidth',2)
hold on
semilogx(Tsw,p_houtinMu(2,:),'rx','MarkerSize',10,'LineWidth',2)

semilogx(Tsw,p_houtinMu(3,:),'ko','MarkerSize',10,'LineWidth',2)
semilogx(Tsw,p_houtinMu(4,:),'ro','MarkerSize',10,'LineWidth',2)

grid on 

% yLim = ylim;
% yLim(1) = 0;
% ylim(yLim)

legend('N Pi-lump','N Pi-lump with unsteady fric.','fMOC','fMOC with unsteady fric.')
ylabel('Pressure, mean (MPa)')
xlabel('Switching period (s)') 


sgtitle('Variation in Pressure and Mean Pressure: High-Pressure Pipeline')
