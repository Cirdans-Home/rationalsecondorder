%% GAUCHI-TYPE INTEGRATOR V1
% Solution of a test wave equation via the Gautschi-type method in time
% and the linear Lagrangian Finite Element in Space. The code uses the
% matlab builtin function to assemble the FEM matrices.

clear; close all hidden; clc;

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

%% Create PDE structure and mesh (for space)
hmax = 0.1;    % Length parameter of the mesh
waveeq = createpde(numberOfPDE);    
waveeq_geom = geometryFromEdges(waveeq, g );
msh = generateMesh(waveeq,'GeometricOrder','quadratic','hmax',hmax);
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

%% Create PDE structure and mesh (for velocity)
waveeqv = createpde(numberOfPDE);    
waveeqv_geom = geometryFromEdges(waveeqv, g );
msh = generateMesh(waveeqv,'GeometricOrder','quadratic','hmax',hmax);
specifyCoefficients(waveeqv,'m', M_COEFF,'d', D_COEFF,'c',...
    C_COEFF,'a', A_COEFF,'f', F_COEFF );

%% Boundary conditions for the first derivative
G_Dprime = @(location,state) uprime(location.x,location.y,state.time);
for ie = 1: NumEdges
    applyBoundaryCondition(waveeqv,'Dirichlet',...
        'edge',ie,'h',1,'r',G_Dprime);
end

%% Time Integration Routine
nt = 1001;
t0 = 0;
t1 = 1;

yzero = @(x,y) utrue(x,y,0).';
yone = @(x,y) uprime(x,y,0).';


options.nt = nt;
options.type = "exact";
options.verbose = 2;
options.ytrue = @(x,y,t) utrue(x,y,t).';
options.yprime = @(x,y,t) uprime(x,y,t).';
[T,Y,options] = gautschint(waveeq,waveeqv,yzero,yone,[t0,t1],options);


