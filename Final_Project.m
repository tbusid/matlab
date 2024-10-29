clear vars
close all;
a = 1;
b = 0;
d = 17;
e = 0;
ode_RHS = @(t,y) [-a*y(1)*y(2)-b*y(1); -d*y(1) - e*y(2)];
T = 10;
N = 1000;
t = 0:T/N:T;
%init_cond = [7000 15000];
init_cond = [10 10 ; 10 5 ; 10 7 ; 5 10 ; 7 10 ; 10 8 ; 8 10 ; 6 10 ; 10 6 ; 10 4 ; 4 10 ; 10 3 ; 3 10];
figure(1);
for i = 1:length(init_cond(:,1))
[tsoln,ysoln] = ode45(ode_RHS,t,init_cond(i,:));
plot(ysoln(:,1),ysoln(:,2),'LineWidth',2,'Color',[0 0.4 0.7]);
hold on;
end
[P,Q] = meshgrid(-2:0.5:10,-2:0.5:10);

U = - a.*P.*Q;% - b.*P;
V = -d.*P;% - e.*Q;

U1 = U./sqrt(U.^2 + V.^2);
V1 = V./sqrt(U.^2 + V.^2);

quiver(P,Q,U1,V1,'Color',[0.7 0 0.4]);

scatter(0,0,100,'or','filled');
scatter(b*e/(a*d),-b/a,100,'or','filled');

set(gca,'FontSize',20);
plot([-5; 20],[0; 0],'k','LineWidth',1);
plot([0; 0],[-5; 20],'k','LineWidth',1);
xlabel('Guerrilla');
ylabel('Conventional');

xlim([-2 10]);
ylim([-2 10]);
print('Final_Project_Graph_13','-dpng');