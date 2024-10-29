close all;
clear;
a = 1;
b = 0;
c = 40;
d = 3;
e = 0;
fplot(@(x)(2*b*x/a+c)^(0.5));

%xlim([-2 10]);
%print('Graphof_Analytical_Solution','-dpng');
ode_RHS = @(t,y) [-a*y(1)*y(2)-b*y(1); -d*y(1) - e*y(2)];
T = 10;
N = 1000;
t = 0:T/N:T;
init_cond = [10 10];
figure(1);
f = @(x)(2*b*x/a+c).^(0.5);
x = 0:0.01:10;
hold on;
for i = 1:length(init_cond(:,1))
[tsoln,ysoln] = ode45(ode_RHS,t,init_cond(i,:));
plot(ysoln(:,1),ysoln(:,2),'LineWidth',2);
end
y = f(ysoln(:,1));
xlim([0 10]);
err = max(abs(ysoln(:,2)-y));
