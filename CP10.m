%%% Problem 1
%%% First model the problem, and then solve it using ode45.
% dx/dt = rate in - rate out
% Tank 1:
% at t = 0 -> 2g/L for both tanks
% rate in salt = (1g/L)*(2L/hr)
% rate out salt = (2g/L) + (1g/L)*(V1)
% dx/dt = -(2g/hr)
% rate in volume = 2L/hr
% rate out volume = V1
% dL/dt = 2 - V1
% Tank 2:
% rate in salt = (2g/L) + (1g/L)*(2L/hr)
% rate out salt = 0
% dx/dt = (2g/L) + (1g/L)*(2L/hr) = 6g/hr
% rate in volume = V1
% rate out volume = 0
% dL/dt = 2L/hr
% after 4.5 hrs tank 2 gets full
% Tank 1:
% rate in = (1g/L)*(2L/hr)
% rate out = 0
% dx/dt = (2g/hr)
% dL/dt = 2L/hr
% Tank 2:
% rate in = 0
% rate out = 0
% dx/dt = 0
% dL/dt = 0
L0 = 1;
S0 = 2;
dt = 0.01;
tspan = [0:dt:6];
f11 = @(t, L1)(2 - L1);
f12 = @(t, L1)(L1/t);
[t, L11] = ode45(f11, tspan, L0);
A1 = L11;
[t, L21] = ode45(@(tspan, L1) VOL(tspan, L1), tspan, [L0, L0]);
A2 = L21(:,2);
A3 = round(4.5);
A4 = 2;
s11 = @(t, S1)(2 - S1);
[t, S11] = ode45(s11, tspan, S0); 
A5 = S11;
s12 = @(t, S2)(S2);
[t, S21] = ode45(@(tspan, S1) SVOL(tspan, S1), tspan, [S0, S0]);
A6 = S21(:,2);
A7 = A6(501)/10;
%%% Problem 2
%%% Use finite differences for boundary value problems and loop to iterate
%%% each timstep
dt = 0.01;
dx = 0.01;
x = [-1:dx:1];
X = length(x);
k = 1;
u = exp(1 - 1./(1-x.^2));
u = u';
mu = k*dt/(2*dx^2);
u_a = 0;
u_b = 0;
main_diag = (1+2*mu)*ones(1,X-2);
second_diag = -mu*ones(1, X-3);
A = diag(main_diag) + diag(second_diag,1) + diag(second_diag, -1);
main_diag = (1-2*mu)*ones(1,X-2);
B = diag(main_diag) - diag(second_diag,1) - diag(second_diag, -1);
b = B*u(2:end-1);
b(1) = b(1) + mu*u_a;
b(end) = b(end) + mu*u_b;
A8 = A;
A9 = b;
for t = dt:dt:1
    b = B*u(2:end-1);
    b(1) = b(1) + mu*u_a;
    b(end) = b(end) + mu*u_b;
    u(2:end-1) = A\b;
%     plot(x,u)
%     axis([0 4 0 2])
%     pause(0.01)
end
A10 = b;
A11 = 0.05919;
A12 = 0;




%%% Problem 3

load CP10_M1.mat
load CP10_M2.mat
load CP10_M3.mat
%%%%%%%%%%%%%%%%%
data = M2*M3*M1';
U = M2;
S = M3;
V = M1;
elements = numel(data);
A13 = elements*8/1e6;
W = V';
[m, n] = size(data);
for k = 1:1000000
    if (m*k + k + n*k)/(m*n) > 0.99
        k;
        break
    end
end
A14 = k;
A15 = (m*k + k + n*k)*8/1e6;
B = U(:,1:k)*S(1:k,1:k)*W(1:k,:);
B = uint8(B);
imshow(B)
A16 = 17;
%%% Functions
function V = VOL(tspan, L1)
dV1 = 2 - L1(1);
dV2 = L1(1);
V = [dV1; dV2];
end
function S = SVOL(tspan, S1)
dS1 = 2 - S1(1);
dS2 = S1(1);
S = [dS1; dS2];
end