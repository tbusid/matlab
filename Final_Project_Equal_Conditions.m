init_cond = [10 10];
b = 5;
a = 2*b/(init_cond(:,2))-0.1;
ode_RHS = @(t,y) [-a*y(1)*y(2); -b*y(1)];
T = 10;
N = 1000;
t = 0:T/N:T;
figure(1);
for i = 1:length(init_cond(:,1))
[tsoln,ysoln] = ode45(ode_RHS,t,init_cond(i,:));
plot(ysoln(:,1),ysoln(:,2),'LineWidth',2,'Color',[0 0.4 0.7]);
hold on;
end
[P,Q] = meshgrid(-2:0.5:10,-2:0.5:10);

U = - a.*P.*Q;
V = -d.*P;

U1 = U./sqrt(U.^2 + V.^2);
V1 = V./sqrt(U.^2 + V.^2);

quiver(P,Q,U1,V1,'Color',[0.7 0 0.4]);

scatter(0,0,100,'or','filled');

set(gca,'FontSize',20);
plot([-5; 20],[0; 0],'k','LineWidth',1);
plot([0; 0],[-5; 20],'k','LineWidth',1);
xlabel('Guerrilla');
ylabel('Conventional');

xlim([-2 10]);
ylim([-2 10]);
print('Final_Project_Graph_15','-dpng');
