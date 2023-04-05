%% GAUCHI-TYPE INTEGRATOR V0.5
% Solution of a test wave equation via the Gautschi-type method in time
% and the linear Lagrangian Finite Element in Space. The code uses the
% matlab builtin function to assemble the FEM matrices.

clear; close all; clc;

%% Problem data
% Test problem with known solution and specific Dirichlet boundary
% conditions on all the boundaries.
f = @(x,y,t) pi^2*(2*t^2 - (x+y).^2).*sin(pi*t.*(x+y));
utrue = @(x,y,t) sin(pi*t.*(x+y));
uprime = @(x,y,t) pi.*(x + y).*cos(pi*t.*(x + y));

%% Geometry
% Circle minus a square domain
gd = [     1     3
    0     4
    0     0
    1     1
    0     1
    0     0
    0     0
    0     0
    0     1
    0     1];
sf = 'C1-SQ1';
ns = [67    83
    49    81
    0    49];
g = decsg(gd,sf,ns);

%% WAVE EQ. COEFFICIENTS
numberOfPDE = 1; % Scalar problem
M_COEFF = 1; % 2nd derivative in time
D_COEFF = 0; % 1st derivative in time
C_COEFF = 1; % Diffusion tensor/coefficient
A_COEFF = 0; % Reaction coefficient
F_COEFF = @(location,state) f(location.x,location.y,state.time); % Source term
G_D = @(location,state) utrue(location.x,location.y,state.time);

%% Create PDE structure and mesh
hmax = 0.05;    % Length parameter of the mesh
waveeq = createpde(numberOfPDE);
waveeq_geom = geometryFromEdges(waveeq, g );
msh = generateMesh(waveeq,'GeometricOrder','linear','hmax',hmax);
[P,E,T] = meshToPet(msh );
figure(1)
subplot(1,2,1)
pdegplot(waveeq,"EdgeLabels","on");
subplot(1,2,2)
pdeplot(P,E,T)
axis tight

specifyCoefficients(waveeq,'m', M_COEFF,'d', D_COEFF,'c',...
    C_COEFF,'a', A_COEFF,'f', F_COEFF );
NumEdges = waveeq_geom.NumEdges;

%% Fix Boundary Conditions (Time dependent)
for ie = 1: NumEdges
    applyBoundaryCondition(waveeq,'Dirichlet','edge',ie,'h',1,'r',G_D);
end

%% Time Domain
t0 = 0;
t1 = 2;
%% Assemble FEM Matrices
% First we assemble the matrices containing the information on the boundary
state.time = t0;
FEMb = assembleFEMatrices(waveeq,'boundary',state);
H = FEMb.H;
[B,Or] = pdenullorth(H);
fem_matrices = assembleFEMatrices(waveeq,'MK',state);
K = fem_matrices.K;
M = B'*fem_matrices.M*B;
Kc = B'*K*B;

%% Initial Conditions
yzero = utrue(P(1,:),P(2,:),0).';
yprime = uprime(P(1,:),P(2,:),0).';

%% Gautschi-Type
mmin = -1;
mmax = 1;
mfunoptions.type = "rational";
mfunoptions.poltype = "PADEEXP";       % Used only if type is rational
mfunoptions.numpoles = 6;

t = [0,2];
h = 0.01;
[T,U] = gautschigen(@(t,y) odefun(t,y,H,Or,B,K,Kc,M,waveeq)...
    ,[t0,t1],yzero,yprime,K,h,mfunoptions);
U = real(U);

%% Plot Solution
figure(3)
pause()
nt = length(T);
for i=1:nt
    subplot(1,3,1)
    pdeplot( waveeq ,'XYData', U( : ,i),...
        'ZData',U(:,i),'ColorMap','jet');
    caxis([mmin,mmax]);
    axis([-1 1 -1 1 mmin mmax]);
    title(sprintf("Time %1.2f",T(i)));
    subplot(1,3,2)
    pdeplot( waveeq ,'XYData', utrue( P(1,:),P(2,:),T(i)).',...
        'ZData',utrue( P(1,:),P(2,:),T(i)).','ColorMap','jet');
    caxis([mmin,mmax]);
    axis([-1 1 -1 1 mmin mmax]);
    subplot(1,3,3)
    pdeplot( waveeq ,'XYData', ....
        abs(U(:,i) - utrue( P(1,:),P(2,:),T(i)).'),'ColorMap','jet');
    pause()
end


%% Odefun function
function y = odefun(t,y,H,Or,B,K,Kc,M,waveeq)
%ODEFUN implements the dynamic for the Gautschi type integrator using the
%MATLAB routines to assemble FEM matrices.

state.time = t;
fem_matrices = assembleFEMatrices(waveeq,'FRG',state);
% Vector with known value of the constraint DoF.
ud = Or*((H*Or\fem_matrices.R));
Fc = B'*(fem_matrices.F + fem_matrices.G - K*ud);
y = B*(M\(-Kc*(B'*y) + Fc))+ud;

end