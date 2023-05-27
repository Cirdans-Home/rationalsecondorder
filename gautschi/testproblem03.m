%% TEST OF GENERAL GAUTSCHI-TYPE INTEGRATOR
% Test a Gautschi type integrator using all matrix function strategies.
% the results are checked against a Matlab ODE solver.

clear; clc; close all;

addpath('../rationalapprox/')

NN = 100;
A = (gallery("toeppen",NN,NN)*gallery("toeppen",NN,NN).');
f = @(t) 0.5*sin(t);
yzero = ones(NN,1);
yprime = zeros(NN,1);

%% Checking order
t = [0,1];
hval = [1e-1,1e-2,1e-3,1e-4];
err = zeros(length(hval),1);
index = 1;
fprintf("  & Gautschi          & ETD2RK            \\\\\n");
fprintf("h      & T (s) & Rel. Err. & T (s) & Rel. Err. \\\\\n");
for h = hval
    %% Solution with Gautschi Solver
    odefun = @(t,y) -A*y + f(t)*ones(NN,1);
    mfunoptions.type = "direct";
    if h >= 1e-2
        mfunoptions.numpoles = 6;
    else
        mfunoptions.numpoles = 2;
    end
    mfunoptions.verbose = false;
    tic;
    [T2,Y2,mfunoptions] = gautschigen(odefun,t,yzero,yprime,A,h,mfunoptions);
    time1 = toc;
    %% Solution with Expint solver
    I = eye(size(A));
    O = zeros(size(A));
    problem.L = [O I;-A O];
    problem.N = @(u,t,problem) [zeros(size(u,1)/2,1);f(t)*ones(size(u,1)/2,1)];
    problem.y0 = [yzero;yprime];
    tic;
    [T,Y] = expglm(problem,t,h,@abnorsett2);
    time2=toc;
    %% Solution with Matlab Solver
    M = [O I;-A O];
    odefun = @(t,y) M*y + [zeros(size(y,1)/2,1);f(t)*ones(size(y,1)/2,1)];

    opts = odeset("AbsTol",1e-12,"RelTol",3e-14);
    [~,Ymat] = ode15s(@(t,y) odefun(t,y),T2,[yzero;yprime],opts);

    %% Computing the errors
    figure(index)
    subplot(2,1,1)
    plot(1:NN,Y(end,1:NN),'r-', ...
        1:NN,Y2(1:NN,end),'b--', ...
        1:NN,Ymat(end,1:NN),'k.','LineWidth',2)
    legend({'ETD2RK','Gautschi','ode15s'},'Location','eastoutside')
    subplot(2,1,2)
    semilogy(1:NN,abs(Y(end,1:NN)-Ymat(end,1:NN)),'r-', ...
        1:NN,abs(Y2(1:NN,end).' - Ymat(end,1:NN)),'b--','LineWidth',2)
    legend({'ETD2RK','Gautschi'},'Location','eastoutside')
    error1 = norm(Ymat(end,1:NN).'- Y2(:,end),2)/norm(Ymat(end,1:NN),2);
    error2 = norm(Ymat(end,1:NN)- Y(end,1:NN),2)/norm(Ymat(end,1:NN),2);

    %% Plotting the Table Lines
    fprintf("%1.1e & %1.2e & %1.2e & %1.2e & %1.2e \\\\\n",h,time1,error1,time2,error2)
    index = index+1;
end