m = 1.0;
load ElCentro.txt;
P=-386.4*m*ElCentro;
dt = 0.02;


Tn = 1.0;
wn = 2*pi/Tn;
k = m*wn^2;
zeta = 0.02;
c=2*zeta*m*wn;

gamma=0.5;
beta=0.25;

[u,v,a]=Newmark(m,c,k,0,0,dt,P,gamma,beta);
%[u,v,a]=InterpolatedLoad(m,c,k,0,0,dt,P);

nSteps=size(P)-1;
t=[0:dt:nSteps*dt];
plot(t,u)
grid

k
m
zeta

U0 = max(abs(u));
F0 = U0*k;

check = @(dU,dR)abs(dR)<1e-12

spring0=ElasticPPSpring(k,250);
[u,v,a,fs]=NonlinearNewmark(m,c,0,0,dt,P,gamma,beta, ...
			    spring0,'Newton',5,check);

Umax = max(abs(u))
Fmax = max(abs(fs))

spring1=ElasticPPSpring(k,125.0);
[u1,v1,a1,fs1]=NonlinearNewmark(m,c,0,0,dt,P,gamma,beta, ...
				spring1,'Newton',5,check);

Umax = max(abs(u1))
Fmax = max(abs(fs1))

spring2=ElasticPPSpring(k,60.0);
[u2,v2,a2,fs2]=NonlinearNewmark(m,c,0,0,dt,P,gamma,beta, ...
				spring2,'Newton',5,check);


Umax = max(abs(u2))
Fmax = max(abs(fs2))

spring3=ElasticPPSpring(k,30.0);
[u3,v3,a3,fs3]=NonlinearNewmark(m,c,0,0,dt,P,gamma,beta, ...
				spring3,'Newton',5,check);


Umax = max(abs(u3))
Fmax = max(abs(fs3))

