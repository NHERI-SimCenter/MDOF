load ElCentro.txt;
P=-386.4*m*ElCentro;
dt = 0.02;

Ts=[.5 1.0 2.0];
for i = 1:3
Tn = Ts(i);
m = 1.0;

wn = 2*pi/Tn;
k = m*wn^2;

wn = sqrt(k/m);
Tn = 2*pi/wn
zeta = 0.0;
c=2*zeta*m*wn;

gamma=0.5;
beta=0.25;

[u,v,a]=Newmark(m,c,k,0,0,dt,P,gamma,beta);
%[u,v,a]=InterpolatedLoad(m,c,k,0,0,dt,P);

nSteps=size(P)-1;
t=[0:dt:nSteps*dt];

%plot(t,u)
%grid

U0 = max(abs(u))

end
