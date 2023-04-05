function [T,Y,options] = gautschint(waveeq,waveeqv,yzero,yone,t,options)
%%GAUCHINT Gauchi-Type Integrator for second order differential equation
% Applies gauchy type integrator to a FEM discretized second order linear
% wave equation implemented in the MATLAB language.
%   INPUT: waveeq Matlab structure for time-dependent FEM problem
%          waveeqprime Matlab structure for time-dependent FEM problem
%          (Velocity BCs)
%          yzero  Initial condition on the state
%          yprime Initial condition on the first derivative
%          t vector of all the time-steps or initial and final times
%          options structure containing the following fields
%           nt      Number of time steps (default 1000)
%           type    Type of application of the matrix functions:
%                   - "exact" Beware: cubic op-cost quadratic cost space
%                   - "ratkrylov" Rational-Krylov with poles determined by
%                       error tollerance
%          If you run the code without inputs it outputs an options
%          variable with all the correct fields inizialized to the default
%          values.

if nargin == 0
    % Default option structure
    options.nt = 1000;
    options.type = "exact";
    options.verbose = true;
    T = [];
    Y = [];
    return
end

if length(t) == 2
    if sum(strcmp(fieldnames(options), 'nt')) == 1
        T = linspace(t(1),t(2),options.nt);
    else
        T = linspace(t(1),t(2),1000);
        options.nt = 1000;
    end
elseif length(t) > 2
    T = t;
else
    t = [0,t];
    gauchint(waveeq,t,options);
end

%% Build the time-independent part of the model
timebuildStart = tic;
state.time = T(1);
FEMb = assembleFEMatrices(waveeq,'boundary',state);
[B,Or] = pdenullorth(FEMb.H);
fem_matrices = assembleFEMatrices(waveeq,'MK',state);
K = fem_matrices.K;
M = fem_matrices.M;
Mc = B'*M*B;
Kc = B'*K*B;
options.timebuildEnd = toc(timebuildStart);
if options.verbose
    fprintf("Time to build FEM matrices:   %1.2e (s)\n",options.timebuildEnd);
    fprintf("Dofs:                         %d\n",size(Kc,1));
    fprintf("nnzs:                         %d\n",nnz(Kc));
end

%% Evaluate initial conditions
P = waveeq.Mesh.Nodes;              % Nodes of the mesh
if options.verbose > 1
    xmmin = min(P(1,:));
    xmmax = max(P(1,:));
    ymmin = min(P(2,:));
    ymmax = max(P(2,:));
end
ndof = size(K,1);
v = zeros(ndof,options.nt);
Y = zeros(ndof,options.nt);
% The elements are Lagrangian therefore the evaluation of initial
% conditions is just the evaluation of the two functions on the nodes
Y(:,1) = yzero(P(1,:),P(2,:));
v(:,1) = yone(P(1,:),P(2,:));

%% March in time

if options.verbose
    wbar = waitbar(0,"Integration underway...");
end
options.rhsassemble = 0;
if sum(strcmp(fieldnames(options), 'ytrue')) == 1
    options.error = zeros(options.nt-1,1);
    options.errorv = zeros(options.nt-1,1);
end

switch upper(options.type)
    case "EXACT"
        h = T(2)-T(1);
        tic;
        psiK = psi(h^2*K); % Matrix function
        options.timemfuneval = toc;
        v(:,1) = sigma(h^2*K)*v(:,1);
        for i=1:options.nt-1
            % Step-length
            if abs((T(i+1)-T(i))-h) > 10*eps
                h = T(i+1)-T(i);
                % Recompute matrix function in case of enough change.
                tic;
                psiK = psi(h^2*K);
                options.timemfuneval = options.timemfuneval+toc;
            end
            % Compute right-hand side
            state.time = T(i);
            tic
            fem_matrices = assembleFEMatrices(waveeq,'FRG',state);
            % Vector with known value of the constraint DoF (Space)
            ud = Or*((FEMb.H*Or\fem_matrices.R));
            Fc = B'*(fem_matrices.F + fem_matrices.G - K*ud);
            % Vector with known value of the constraint DoF (Velocity)
            state.time = T(i) + 0.5*h;
            fem_matrices = assembleFEMatrices(waveeqv,'R',state);
            vd = Or*((FEMb.H*Or\fem_matrices.R));
            options.rhsassemble = options.rhsassemble + toc;
            % half-step on velocity
            ytemp = B*(Mc\(-Kc*(B'*Y(:,i)) + Fc))+ud;
            vtemp = B*(B'*(v(:,i) + 0.5*h*psiK*ytemp))+vd;
            % Compute new right-hand side
            tic;
            state.time = T(i+1);
            fem_matrices = assembleFEMatrices(waveeq,'FRG',state);
            % Vector with known value of the constraint DoF.
            ud = Or*((FEMb.H*Or\fem_matrices.R));
            Fc = B'*(fem_matrices.F + fem_matrices.G - K*ud);
            fem_matrices = assembleFEMatrices(waveeqv,'R',state);
            vd = Or*((FEMb.H*Or\fem_matrices.R));
            options.rhsassemble = options.rhsassemble + toc;
            % new-step on space
            Y(:,i+1) = B*(B'*(Y(:,i) + h*vtemp))+ud;
            % full-step on velocity
            ytemp = B*(Mc\(-Kc*(B'*Y(:,i+1)) + Fc))+ud;
            v(:,i+1) = B*(B'*(vtemp + 0.5*h*psiK*ytemp)) + vd;
            %v(:,i+1) = vtemp + 0.5*h*psiK*ytemp;
            if options.verbose
                wbar = waitbar(((i+1)/(options.nt+1)),wbar,...
                    "Integration underway...");
            end
            if options.verbose > 1
                figure(42)
                if sum(strcmp(fieldnames(options), 'ytrue')) == 1
                    subplot(2,2,1)
                else
                    subplot(1,2,1)
                end
                pdeplot( waveeq ,'XYData', Y( : ,i+1),...
                    'ZData',Y(:,i+1),'ColorMap','jet');
                title("Numerical solution")
                xlim([xmmin xmmax]); 
                ylim([ymmin ymmax]);
                if sum(strcmp(fieldnames(options), 'ytrue')) == 1
                    subplot(2,2,2)
                else
                    subplot(1,2,2)
                end
                pdeplot( waveeq ,'XYData', v( : ,i+1),...
                    'ZData',v(:,i+1),'ColorMap','jet');
                title("Numerical velocity")
                if sum(strcmp(fieldnames(options), 'ytrue')) == 1
                    subplot(2,2,3)
                    pdeplot( waveeq ,...
                        'XYData', options.ytrue( P(1,:),P(2,:),T(i+1)),...
                        'ZData',options.ytrue( P(1,:),P(2,:),T(i+1)),...
                        'ColorMap','jet');
                    title("True solution")
                    subplot(2,2,4)
                    pdeplot( waveeq ,...
                        'XYData', options.yprime( P(1,:),P(2,:),T(i+1)),...
                        'ZData',options.yprime( P(1,:),P(2,:),T(i+1)),...
                        'ColorMap','jet');
                    title("True velocity")
                    xlim([xmmin xmmax]); 
                    ylim([ymmin ymmax]);
                    figure(43)
                    subplot(1,2,1)
                    pdeplot( waveeq ,...
                        'XYData',...
                        abs(Y( : ,i+1) - options.ytrue( P(1,:),P(2,:),T(i+1))),...
                        'ZData',...
                        abs(Y( : ,i+1) - options.ytrue( P(1,:),P(2,:),T(i+1))),...
                        'ColorMap','jet');
                    title('Absolute Error (space)')
                    xlim([xmmin xmmax]);
                    ylim([ymmin ymmax]);
                    subplot(1,2,2)
                    pdeplot( waveeq ,...
                        'XYData',...
                        abs(v( : ,i+1) - options.yprime( P(1,:),P(2,:),T(i+1))),...
                        'ZData',...
                        abs(v( : ,i+1) - options.yprime( P(1,:),P(2,:),T(i+1))),...
                        'ColorMap','jet');
                    title('Absolute Error (velocity)')
                    xlim([xmmin xmmax]);
                    ylim([ymmin ymmax]);
                end
                pause(1/25)
            end
            if sum(strcmp(fieldnames(options), 'ytrue')) == 1
                % Compute error
                options.error(i) = norm(Y( : ,i+1) - ...
                    options.ytrue( P(1,:),P(2,:),T(i+1)),2);
                options.errorv(i) = norm(v(:,i+1) - ...
                    options.yprime(P(1,:),P(2,:),T(i+1)),2);
                fprintf("At t = %e ||y_n - y(t_n)|| = %e " + ...
                    "||v_n - y'(t_n)|| = %e\n",T(i+1),options.error(i),...
                    options.errorv(i));
                if options.error(i) > 1
                    figure(404)
                    subplot(1,3,1)
                    pdeplot( waveeq ,'XYData', Y( : ,i+1),...
                        'ZData',Y(:,i+1),'ColorMap','jet');
                    title("Numerical solution")
                    subplot(1,3,2)
                    pdeplot( waveeq ,'XYData', options.ytrue( P(1,:),P(2,:),T(i+1)),...
                        'ZData',options.ytrue( P(1,:),P(2,:),T(i+1)),'ColorMap','jet');
                    title("Numerical solution")
                    subplot(1,3,3)
                    pdeplot( waveeq ,'XYData', abs(Y( : ,i+1)-options.ytrue( P(1,:),P(2,:),T(i+1))),...
                        'ZData',abs(Y(:,i+1)-options.ytrue( P(1,:),P(2,:),T(i+1))),'ColorMap','jet');
                    title("Error")
                    pause()
                end
            end
        end
    case "RATKRYLOV"
    otherwise
        errorr("This is an uknown method!")
end

if options.verbose
    fprintf("Dense matrix value evalution: %1.2e (s)\n",options.timemfuneval);
    fprintf("Time to assemble rhs vectors: %1.2e (s)\n",options.rhsassemble);
    fprintf("\n\n")
    close(wbar);
end

end

%% Auxiliary routines for the computation of the matrix function
% These two routines computes the dense matrix function, should be used
% only for small problem and/or to perform against the matrix-vector
% product problems.
function S = sigma(A)
%%SIGMA computes S = sqrt(A)^{-1} sin(sqrt(A))
sqA = sqrtm(full(A));
S = real(sqA\funm(sqA,@sin));
end

function P = psi(A)
%%PSI computes P = ( (sqrt(A)/2)^{-1} sin(sqrt(A)/2) )^2
sqA2 = sqrtm(full(A))/2;
P = sqA2\funm(sqA2,@sin);
P = real(P*P);
end