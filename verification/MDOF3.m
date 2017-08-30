clear all;

load ElCentro.txt
g = 386.4;
ug = g*ElCentro;
t=[0:.02:1559*.02]';
dt = 0.02;

m=100/g;
k=31.54;
%m =10

numStory = 5;

alpha = 0.01;
beta = 0.001;

zeta = 0.05;

M=eye(numStory)*m; 
M(numStory,numStory) = alpha*m;

K=eye(numStory)*2*k; 
for i=1:numStory-2
    K(i+1,i)=-k;
    K(i,i+1)=-k;
end

K(numStory-1,numStory-1)=k+beta*k;
K(numStory,numStory)=beta*k;
K(numStory-1,numStory)=-beta*k;
K(numStory,numStory-1)=-beta*k;

M*386.4
K

[eigVectors, eigValues]=eig(K,M);

eigValues=diag(eigValues);
natFrequencies=sqrt(eigValues);
eigVectors;

periods=ones(numStory,1)*2*pi;
periods=periods./natFrequencies

numDamp = 3;
fsMax=zeros(numStory,numDamp+1);
fdMax=zeros(numStory,numDamp+1);

w1=natFrequencies(1);
wb=natFrequencies(3);

C=zeros(numStory,numStory);
for i=1:numStory
  C=C+2*zeta*natFrequencies(i)*eigVectors(:,i)*eigVectors(:,i)';
end
C = M*C*M;

P=-M*ones(numStory,1)*ug';
    
gamma = 0.5;
beta=0.25;
[u,v,a] = Newmark(M,C,K,zeros(numStory,1),zeros(numStory,1),dt,P,gamma,beta);

max(abs(u(numStory,:)))
