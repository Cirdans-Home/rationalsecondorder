%% TEST OF GENERAL GAUTSCHI-TYPE INTEGRATOR
% Test a Gautschi type integrator using dense matrix function computations.
% the results are checked against a Matlab ODE solver.

clear; clc; close all;

A = [2 -1 0; -1 2 0; 0 -1 2];
f = @(t) 0.5*sin(t);
yzero = rand(3,1);
yprime = zeros(3,1);

%% Solution with Gautschi Solver
odefun = @(t,y) -A*y + f(t)*ones(3,1);
t = [0,2];
h = 0.000001;
mfunoptions.type = "rational";
mfunoptions.poltype = "PADEEXP";       % Used only if type is rational
mfunoptions.numpoles = 3;
[T2,Y2,mfunoptions] = gautschigen(odefun,[0,2],yzero,yprime,A,h,mfunoptions);

%% Solution with Matlab solver
I = eye(size(A));
O = zeros(size(A));
M = [O I;-A O];
odefun = @(t,y) M*y + [zeros(3,1);f(t)*ones(3,1)];

opts = odeset("AbsTol",1e-10,"RelTol",1e-13);
[T,Y] = ode15s(@(t,y) odefun(t,y),T2,[yzero;yprime],opts);

figure(1)
plot3(Y(:,1),Y(:,2),Y(:,3),'b-.',...
    Y2(1,:),Y2(2,:),Y2(3,:),'r--','LineWidth',2)

figure(2)
subplot(1,3,1)
plot(T,Y(:,1),'b-.',T2,Y2(1,:),'r--','LineWidth',2);
legend('ode45','Gautschi')
subplot(1,3,2)
plot(T,Y(:,2),'b-.',T2,Y2(2,:),'r--','LineWidth',2);
legend('ode45','Gautschi')
subplot(1,3,3)
plot(T,Y(:,3),'b-.',T2,Y2(3,:),'r--','LineWidth',2);
legend('ode45','Gautschi')

for i=1:length(T)
    err(i) = norm(Y(i,1:3).' -  Y2(:,i),2);
end
figure(3)
semilogy(T,err,'k','LineWidth',2)