%% GAUCHI-TYPE INTEGRATOR FOR THE WAVE EQUATION
% Solution of a test wave equation via the Gautschi-type method in time
% and the linear Lagrangian Finite Element in Space. The code uses the
% matlab builtin function to assemble the FEM matrices.

clear; close all; clc;

addpath('../rationalapprox/')

%% Problem data
% Test problem with known solution and specific Dirichlet boundary
% conditions on all the boundaries.
f = @(x,y,t) 0*x + 0*y + 0*t;
uzero = @(x,y,t) 0.8*exp(-(x+0.3).^2/0.06 -(y+0.3).^2/0.06) + 0*t;
uprime = @(x,y,t) 0*x + 0*y + 0*t;

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
G_D = @(location,state) 0*uzero(location.x,location.y,state.time);

%% Create PDE structure and mesh
hmax = 0.03;    % Length parameter of the mesh
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
t1 = 1;
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
yzero = uzero(P(1,:),P(2,:),0).';
yprime = uprime(P(1,:),P(2,:),0).';

%% Gautschi-Type
mmin = -1;
mmax = 1;
[~,~,mfunoptions] = gautschigen();     % Initialize empty structure
mfunoptions.type = "expsum";           % rational, direct, expsum
mfunoptions.poltype = "PADEEXP"; % Used only if type is rational PADEEXP or PADEHYPERGEOM
mfunoptions.numpoles = 6;
mfunoptions.expsumterms = 12;

t = [t0,t1];
h = 1e-2;
[T,U] = gautschigen(@(t,y) odefunrid(t,y,H,Or,B,K,Kc,M,waveeq)...
    ,t,B'*yzero,B'*yprime,Kc,h,mfunoptions);
U = real(U);

%% Check with MATLAB PDE Solver
doplot = false;
if doplot
    Tpde = T;
else
    Tpde = linspace(t0,t1,50);
end
i0 = @(location) uzero(location.x,location.y,0);
i1 = @(location) uprime(location.x,location.y,0);
setInitialConditions(waveeq,i0,i1);
result = solvepde(waveeq,Tpde);

%% Plot Solution
fem_matrices_R = assembleFEMatrices(waveeq,'R',state);
ud = Or*((H*Or\fem_matrices_R.R));
nt = length(T);
if doplot
    figure(3)
    pause()
    err = zeros(nt,1);
    for i=1:nt
        subplot(1,3,1)
        pdeplot( waveeq ,'XYData', B*U( : ,i)+ud,...
            'ZData',B*U(:,i)+ud,'ColorMap','jet');
        caxis([mmin,mmax]);
        axis([-1 1 -1 1 mmin mmax]);
        title(sprintf("Time %1.2f",T(i)));
        subplot(1,3,2)
        pdeplot( waveeq ,'XYData', result.NodalSolution(:,i),...
            'ZData',result.NodalSolution(:,i),'ColorMap','jet');
        caxis([mmin,mmax]);
        axis([-1 1 -1 1 mmin mmax]);
        subplot(1,3,3)
        pdeplot( waveeq ,'XYData', ....
            log10(abs(B*U(:,i)+ud - result.NodalSolution(:,i))),'ColorMap','jet');
        err(i) = norm(B*U(:,i)+ud - result.NodalSolution(:,i))./norm(result.NodalSolution(:,i));
        pause(1/25)
    end
end

%% Error
figure(4)
i = nt;
pdeplot( waveeq ,'XYData', ....
    log10((abs(B*U(:,i) + ud - result.NodalSolution(:,end)))),...
    'ColorMap','jet','Mesh','on');
err_final = norm(B*U(:,i) + ud - result.NodalSolution(:,end));
axis([-1 1 -1 1 mmin mmax]);
xlabel('x');
ylabel('y');
axis square
set(gcf,'Color','white');
export_fig(sprintf('wave_error_%s_%s_%d.eps',...
   mfunoptions.type,mfunoptions.poltype,mfunoptions.numpoles));


%% Odefun function
function y = odefunrid(t,y,H,Or,B,K,Kc,M,waveeq)
%ODEFUN implements the dynamic for the Gautschi type integrator using the
%MATLAB routines to assemble FEM matrices.

state.time = t;
fem_matrices = assembleFEMatrices(waveeq,'FRG',state);
% Vector with known value of the constraint DoF.
ud = Or*((H*Or\fem_matrices.R));
Fc = B'*(fem_matrices.F + fem_matrices.G - K*ud);
y = M\(-Kc*y + Fc);

end