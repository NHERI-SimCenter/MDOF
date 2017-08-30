function [U,V,A,FS] = NonlinearNewmark(M,C,U0,V0,dt,Pt,gamma,beta,...
  structure, algorithm, maxIter, convergenceCheck)
% function to compute response of nonlinear system using Newmark algorithm
%
% Inputs: M  - mass
%         C  - damping coeficient
%         U0 - initial displacement
%         V0 - initial velocity
%         dt - time step, deltaT
%         P  - discretized forcing function, dt time steps P(1) = P0
%         gamma - newmark velocity coeficient
%         beta  - Newmark acceleration coeficient
%         structure -
%         algorithm 'Initial','Modified'or 'Newton'
%         maxIter - max number of iterations per step
%         convergenceCheck - function to evaluate convergence given dU and dR
%
% Outputs: U - displacement vector, discretized dt, U(1) = U0
%          V - velocity vector,     discretized dt, V(1) = V0
%          A - acceleration vector, discretized dt, A(1) = A0
%          Fs - force in structure
%
% written: fmk 09/2016

% Allocation of Returns Values
nSteps=length(Pt);
U=zeros(1,nSteps);
V=zeros(1,nSteps);
A=zeros(1,nSteps);
FS=zeros(1,nSteps);

% form tangent if initial stiffness iterations
[K Fs] = structure.setTrialDisplacement(U0);
if strcmp(algorithm,'Initial')
    Khat = M/(beta*dt^2) + C*gamma/(beta*dt) + K;
end

% set t0 values
U(1)=U0;
V(1)=V0;
A(1)=(Pt(1)-C*V(1)-Fs)/M;

[K Fs] = structure.setTrialDisplacement(U(1));

% loop over all time steps computing response
for i=1:nSteps-1
    
    % form tangent if Modified Newton
    if strcmp(algorithm,'Modified')
        Khat = M/(beta*dt^2) + C*gamma/(beta*dt) + K;
    end
    
    % Predict
    U(i+1)=U(i);
    V(i+1)=(1-gamma/beta)*V(i)+dt*(1.0-gamma/(2*beta))*A(i);
    A(i+1)= (-1/(beta*dt))*V(i)+(1.0-1.0/(2*beta))*A(i);
    dR=Pt(i+1)-M*A(i+1)-C*V(i+1)-Fs; 
    
    % iterate until convergence
    for j=1:maxIter
        
        % Correct
        
        % form tangent if Newton-Raphson
        if strcmp(algorithm,'Newton')
            Khat = M/(beta*dt^2) + C*gamma/(beta*dt) + K;
        end
        
        dU=dR/Khat;
        U(i+1)=U(i+1)+dU;
        V(i+1)=V(i+1)+gamma/(beta*dt)*dU;
        A(i+1)=A(i+1)+1.0/(beta*dt^2)*dU;
        [K Fs] = structure.setTrialDisplacement(U(i+1));
        dR=Pt(i+1)-M*A(i+1)-C*V(i+1)-Fs;
        ok = convergenceCheck(dU,dR);
        if (ok)
            break;
        end;
    end

    % save forces
    FS(i+1) = Fs;

    if (ok)
        structure.commitState();
    else % cannot converge
        break;
    end
    
end


