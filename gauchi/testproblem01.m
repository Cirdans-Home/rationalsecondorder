%% TEST OF GENERAL GAUTSCHI-TYPE INTEGRATOR
% Test a Gautschi type integrator using all matrix function strategies.
% the results are checked against a Matlab ODE solver.

clear; clc; close all;

addpath('../rationalapprox/')

A = (gallery("toeppen",20,20)*gallery("toeppen",20,20).');
f = @(t) 0.5*sin(t);
yzero = ones(20,1);
yprime = zeros(20,1);
lmax = eigs(A,1,"largestabs");

%% Checking order
t = [0,1];
hval = [1e-1,1e-2,1e-3,1e-4];
err = zeros(length(hval),1);
index = 1;
for h = hval
    %% Solution with Gautschi Solver
    odefun = @(t,y) -A*y + f(t)*ones(20,1);
    mfunoptions.type = "direct";
    [T2,Y2,mfunoptions] = gautschigen(odefun,t,yzero,yprime,A,h,mfunoptions);

    %% Solution with Matlab solver
    I = eye(size(A));
    O = zeros(size(A));
    M = [O I;-A O];
    odefun = @(t,y) M*y + [zeros(20,1);f(t)*ones(20,1)];

    opts = odeset("AbsTol",1e-12,"RelTol",3e-14);
    [T,Y] = ode15s(@(t,y) odefun(t,y),T2,[yzero;yprime],opts);

    err(index) = norm(Y(end,1:20).'- Y2(:,end),2)/norm(Y(end,1:20),2);
    index = index + 1;
end
%% Error Plot
hval = hval.';
figure(1)
loglog(hval,err,'ko-',hval,hval.^2,'k.--','LineWidth',2);
legend({'Gautschi','Order 2'},'Location','northeast')
xlabel('h')
ylabel('Relative Error')
axis tight
grid on
set(gca,'xdir','reverse')

%% USING MATRIX FUNCTION ROUTINES
t = [0,1];
hval = [1e-1,1e-2,1e-3,1e-4];
errrat1 = zeros(length(hval),1);
numpoles = zeros(length(hval),1);
index = 1;
for h = hval
    %% Solution with Gautschi Solver
    odefun = @(t,y) -A*y + f(t)*ones(20,1);
    mfunoptions.type = "rational";
    mfunoptions.poltype = "padehypergeom"; % padeexp or padehypergeom
    if strcmpi(mfunoptions.poltype,"PADEEXP")
        mfunoptions.numpoles = numpadpoles(h^2,h^2*lmax);
    elseif strcmpi(mfunoptions.poltype,"PADEHYPERGEOM")
        mfunoptions.numpoles = numlegpoles(h^2,h^2*lmax);
    end
    [T2,Y2,mfunoptions] = gautschigen(odefun,t,yzero,yprime,A,h,mfunoptions);
    numpoles(index) = mfunoptions.numpoles;

    %% Solution with Matlab solver
    I = eye(size(A));
    O = zeros(size(A));
    M = [O I;-A O];
    odefun = @(t,y) M*y + [zeros(20,1);f(t)*ones(20,1)];

    opts = odeset("AbsTol",1e-12,"RelTol",3e-14);
    [T,Y] = ode15s(@(t,y) odefun(t,y),T2,[yzero;yprime],opts);

    errrat1(index) = norm(Y(end,1:20).'- Y2(:,end),2)/norm(Y(end,1:20),2);
    index = index + 1;
end
%% Error Plot
hval = hval.';
figure(2)
loglog(hval,errrat1,'ko-',hval,hval.^2,'k.--','LineWidth',2);
for i = 1:length(hval)
    str = sprintf('  %d Poles',numpoles(i));
    text(hval(i),errrat1(i),str);
end
legend({'Gautschi','Order 2'},'Location','northeast')
xlabel('h')
ylabel('Relative Error')
axis tight
grid on
set(gca,'xdir','reverse')

%% USING EXPONENTIAL SUMS
t = [0,1];
hval = [1e-1,1e-2,1e-3,1e-4];
errrat2 = zeros(length(hval),1);
numpoles = zeros(length(hval),1);
index = 1;
for h = hval
    %% Solution with Gautschi Solver
    odefun = @(t,y) -A*y + f(t)*ones(20,1);
    mfunoptions.type = "expsum";
    mfunoptions.poltype = "padehypergeom"; % padeexp or padehypergeom
    mfunoptions.numpoles = 10;
    mfunoptions.expsumterms = 10;
    [T2,Y2,mfunoptions] = gautschigen(odefun,t,yzero,yprime,A,h,mfunoptions);
    numpoles(index) = mfunoptions.numpoles;

    %% Solution with Matlab solver
    I = eye(size(A));
    O = zeros(size(A));
    M = [O I;-A O];
    odefun = @(t,y) M*y + [zeros(20,1);f(t)*ones(20,1)];

    opts = odeset("AbsTol",1e-12,"RelTol",3e-14);
    [T,Y] = ode15s(@(t,y) odefun(t,y),T2,[yzero;yprime],opts);

    errrat2(index) = norm(Y(end,1:20).'- Y2(:,end),2)/norm(Y(end,1:20),2);
    index = index + 1;
end
%% Error Plot
hval = hval.';
figure(3)
loglog(hval,errrat2,'ko-',hval,hval.^2,'k.--','LineWidth',2);
legend({'Gautschi','Order 2'},'Location','northeast')
xlabel('h')
ylabel('Relative Error')
axis tight
grid on
set(gca,'xdir','reverse')

%% AUXILIARY FUNCTIONS
function n = numpadpoles(h,lmax)
%NUMPADEPOLES estimate of the number of poles needed to reach a given
%accuracy for the approximation of the approximation of the Sinc function
%based on the Padé approximation of the exponential.
z = linspace(0,lmax,ceil(lmax));

for n = 1:15
    tolz = 2*(pi*2^(-4*n)*exp(-z.^2./(8*n+4)).*z.^2.*abs(cos(z./(8*n+4))))./...
        ((2*n+1).*gamma(n+0.5).^2);
    tol = norm(tolz,"inf");
    fprintf("n = %d tol = %e h = %e \n",n,tol,h);
    if tol <= h
        return;
    end
end
end

function n = numlegpoles(h,lmax,z)
%NUMLEGPOLES estimate of the number of poles needed to reach a given
%accuracy for the approximation of the approximation of the Sinc function
%based on the Padé approximation of the Hypergeometric function.

if ~exist("z","var")
    z = linspace(0,lmax,ceil(lmax));
end

for n = 1:15
    tolz = (2*pi*4^(-n-1)*exp(- z.^2/(2*n+2) + (2*n+1)*log(z)))./...
        (gamma(n + 3/2)^2);
    tol = norm(tolz,"inf");
    fprintf("n = %d tol = %e h = %e \n",n,tol,h);
    if tol <= h
        return;
    end
end

end