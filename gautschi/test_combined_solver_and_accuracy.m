%% GAUCHI-TYPE INTEGRATOR FOR THE WAVE EQUATION - Test Time
% Solution of a test wave equation via the Gautschi-type method in time
% and the linear Lagrangian Finite Element in Space. The code uses the
% matlab builtin function to assemble the FEM matrices.

clear; close all; clc;

addpath('../rationalapprox/')
addpath('../chebfun/')
addpath('../rktoolbox/')
addpath('../mkl')

fid = fopen('times_and_accuracy.txt','a+');

fprintf(fid,'\n\nExperiment launched on %s\n\n',datetime);

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

for hmax = [0.1 0.05 0.005 0.0005]  % Length parameter of the mesh %
    %% Create PDE structure and mesh
    waveeq = createpde(numberOfPDE);
    waveeq_geom = geometryFromEdges(waveeq, g );
    msh = generateMesh(waveeq,'GeometricOrder','linear','hmax',hmax);
    [P,E,T] = meshToPet(msh);

    specifyCoefficients(waveeq,'m', M_COEFF,'d', D_COEFF,'c',...
        C_COEFF,'a', A_COEFF,'f', F_COEFF );
    NumEdges = waveeq_geom.NumEdges;

    %% Fix Boundary Conditions (Time dependent)
    for ie = 1: NumEdges
        applyBoundaryCondition(waveeq,'Dirichlet','edge',ie,'h',1,'r',G_D);
    end

    %% Time Domain
    t0 = 0;
    t1 = 0.5;
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
    fem_matrices_R = assembleFEMatrices(waveeq,'R',state);
    ud = Or*((H*Or\fem_matrices_R.R));

    waveeq.SolverOptions.AbsoluteTolerance = 1e-12;
    waveeq.SolverOptions.RelativeTolerance = 1e-12;

    %% Initial Conditions
    yzero = uzero(P(1,:),P(2,:),0).';
    yprime = uprime(P(1,:),P(2,:),0).';

    %% Gautschi-Type
    mmin = -1;
    mmax = 1;
    [~,~,mfunoptions] = gautschigen();     % Initialize empty structure
    mfunoptions.type = "rational";           % rational, direct, expsum
    mfunoptions.poltype = "PADEHYPERGEOM"; % Used only if type is rational PADEEXP or PADEHYPERGEOM
    mfunoptions.expsumterms = 12;

    %% Check with MATLAB PDE solver

    t = [t0,t1];
    lmax = eigs(Kc,1,"largestabs");
    z = linspace(0,lmax,5000);
    h = 1e-1*hmax;
    fprintf(fid,'%1.4f & %1.2e ',hmax,h);
    
    for solvename = ["gmres","pgmres","bicgstab","pbicgstab"]

        mfunoptions.solvetype = solvename;

        mfunoptions.numpoles = numlegpoles(h,h^2*lmax,z);
        tic;
        [T,U] = gautschigen(@(t,y) odefunrid(t,y,H,Or,B,K,Kc,M,waveeq)...
            ,t,B'*yzero,B'*yprime,Kc,h,mfunoptions);
        U = real(U);
        time_gautschi = toc;

        i0 = @(location) uzero(location.x,location.y,0);
        i1 = @(location) uprime(location.x,location.y,0);
        setInitialConditions(waveeq,i0,i1);
        result = solvepde(waveeq,T);

        fprintf(fid,'& %1.2e (%d) & %1.2e ',...
            time_gautschi,mfunoptions.numpoles,...
            max(max(((abs(B*U(:,end) + ud - result.NodalSolution(:,end)))))));
    end
    %% Polynomial solution
    tic;
    [Tp,Up,krylsize] = polgautschi(@(t,y) odefunrid(t,y,H,Or,B,K,Kc,M,waveeq)...
        ,t,B'*yzero,B'*yprime,Kc,h);
    Up = real(Up);
    time_pol = toc;

    i0 = @(location) uzero(location.x,location.y,0);
    i1 = @(location) uprime(location.x,location.y,0);
    setInitialConditions(waveeq,i0,i1);
    result = solvepde(waveeq,Tp);

    fprintf(fid,'& %1.2e [%d] & %1.2e',...
        time_pol,round(mean(krylsize)),...
        max(max(((abs(B*U(:,end) + ud - result.NodalSolution(:,end)))))));

    fprintf(fid,'\\\\\n');
end

fclose(fid);


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

%% Adaptive pole choice
function n = numlegpoles(h,lmax,z)
%NUMLEGPOLES estimate of the number of poles needed to reach a given
%accuracy for the approximation of the approximation of the Sinc function
%based on the Padé approximation of the Hypergeometric function.

if ~exist("z","var")
    z = linspace(0,lmax,ceil(lmax));
end

for n = 1:15
    % Non-symmetric version
    %tolz = (2^(n+1)*factorial(n)/factorial(2*n+1))^2.*z.^(2*n+1);
    % Symmetric version
    tolz = (n+1)/(4*n+6)*(factorial(n)/factorial(2*n+1))^2.*z.^(2*n+2);
    tol = norm(tolz,"inf");
    % fprintf("n = %d tol = %e h = %e \n",n,tol,h);
    if tol <= h
        return;
    end
end

end
