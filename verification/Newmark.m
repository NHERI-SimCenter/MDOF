function [U,V,A] = Newmark(M,C,K,U0,V0,dt,Pt,gamma,beta)
% function to compute response of linear system using Newmark algorithm
%
% Inputs: M  - mass
%         C  - damping coeficient
%         K  - linear stiffness
%         U0 - initial displacement 
%         V0 - initial velocity
%         dt - time step, deltaT
%         P  - discretized forcing function, dt time steps P(1) = P0
%         gamma - newmark velocity coeficient
%         beta  - Newmark acceleration coeficient
%
% Outputs: U - displacements, discretized dt, U(1,:) = U0
%          V - velocities,  discretized dt, V(1,:) = V0
%          A - accelerations, discretized dt, A(1,:) = A0
%
% written: fmk 09/2016

% Allocation of Returns Values
numDOF = size(K,1);

if (size(Pt,2) == numDOF)
    Pt = Pt';  % in case in SOOF model Pt passed by col!
end;

nSteps = size(Pt,2);
U=zeros(numDOF,nSteps);
V=zeros(numDOF,nSteps);
A=zeros(numDOF,nSteps);

% Set Initial values at t0
warning('off','all');
A0 = M\(Pt(:,1)-C*V0-K*U0); 
A0(isnan(A0))=0;% needed for MDOF
A0(isinf(A0))=0;% needed for MDOF
warning('on','all');

size(U)
U(:,1)=U0;
V(:,1)=V0;
A(:,1)=A0;

% Set constants
Khat = 1/(beta*dt^2)*M + gamma/(beta*dt)*C  + K;

% loop over all time steps computing response
for i=1:nSteps-1
    
    % Predict
    U(:,i+1)=U(:,i);
    V(:,i+1)=(1-gamma/beta)*V(:,i)+dt*(1.0-gamma/(2*beta))*A(:,i);
    A(:,i+1)= (-1/(beta*dt))*V(:,i)+(1.0-1.0/(2*beta))*A(:,i);
    
    % Correct
    dR=Pt(:,i+1)-M*A(:,i+1)-C*V(:,i+1)-K*U(:,i+1);
    dU=Khat\dR;
    
    U(:,i+1)=U(:,i+1)+dU;
    V(:,i+1)=V(:,i+1)+gamma/(beta*dt)*dU;
    A(:,i+1)=A(:,i+1)+1.0/(beta*dt^2)*dU;
end


