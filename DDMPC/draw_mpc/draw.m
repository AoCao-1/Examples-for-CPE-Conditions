load('y_cl_1.mat')
load('y_cl_2.mat')
load('y_cl_3.mat')
load('u_cl_1.mat')
load('u_cl_2.mat')
load('u_cl_3.mat')


   

T = 100;
subplot(2,2,1)
plot(0:T-1,y_cl_1(1,1:end),'--','LineWidth',1,'color','green')
hold on
plot(0:T-1,y_cl_2(1,1:end),'--','LineWidth',1,'color','black')
hold on
plot(0:T-1,y_cl_3(1,1:end),'-','LineWidth',1,'color',[0.3,0.75,0.93])
hold on
plot(0:T-1,0*ones(1,T),'-','LineWidth',1,'color','red')
xlabel('Time steps','Fontname','Times New Roman');
ylabel('States','Fontname','Times New Roman');
legend('$x^1$(CCPE)','$x^1$(MCPE)','$x^1$(HCPE)','$x^*$','Interpreter','latex')
ylim([-2 2])

subplot(2,2,2)
plot(0:T-1,y_cl_1(2,1:end),'--','LineWidth',1,'color','green')
hold on
plot(0:T-1,y_cl_2(2,1:end),'--','LineWidth',1,'color','black')
hold on
plot(0:T-1,y_cl_3(2,1:end),'-','LineWidth',1,'color',[0.3,0.75,0.93])
hold on
plot(0:T-1,0*ones(1,T),'-','LineWidth',1,'color','red')
xlabel('Time steps','Fontname','Times New Roman');
ylabel('States','Fontname','Times New Roman');
legend('$x^2$(CCPE)','$x^2$(MCPE)','$x^2$(HCPE)','$x^*$','Interpreter','latex')


subplot(2,2,3)
plot(0:T-1,y_cl_1(3,1:end),'--','LineWidth',1,'color','green')
hold on
plot(0:T-1,y_cl_2(3,1:end),'--','LineWidth',1,'color','black')
hold on
plot(0:T-1,y_cl_3(3,1:end),'-','LineWidth',1,'color',[0.3,0.75,0.93])
hold on
plot(0:T-1,0*ones(1,T),'-','LineWidth',1,'color','red')
xlabel('Time steps','Fontname','Times New Roman');
ylabel('States','Fontname','Times New Roman');
legend('$x^3$(CCPE)','$x^3$(MCPE)','$x^3$(HCPE)','$x^*$','Interpreter','latex')


subplot(2,2,4)
plot(0:T-1,y_cl_1(4,1:end),'--','LineWidth',1,'color','green')
hold on
plot(0:T-1,y_cl_2(4,1:end),'--','LineWidth',1,'color','black')
hold on
plot(0:T-1,y_cl_3(4,1:end),'-','LineWidth',1,'color',[0.3,0.75,0.93])
hold on
plot(0:T-1,0*ones(1,T),'-','LineWidth',1,'color','red')
xlabel('Time steps','Fontname','Times New Roman');
ylabel('States','Fontname','Times New Roman');
legend('$x^4$(CCPE)','$x^4$(MCPE)','$x^4$(HCPE)','$x^*$','Interpreter','latex')