i = 4;
j = 1;


figure(1)

subplot(3,1,1)
plot(out(i,j).t,out(i,j).p_hin)
hold on

subplot(3,1,2)
plot(out(i,j).t,out(i,j).p_hout)
hold on

subplot(3,1,3)
plot(out(i,j).t,out(i,j).p_lin)
hold on

figure(2)
plot(out(i,j).t,out(i,j).q_p)
hold on

title('Pump Flow vs. Time')
xlabel('Time (s)')
ylabel('Flow rate (m^3/s)')

legend('Steady, N Pi-lump with unsteady fric.','Steady, fMOC with unsteady fric.',...
    'Unsteady, N Pi-lump with unsteady fric.','Unsteady, fMOC with unsteady fric.')