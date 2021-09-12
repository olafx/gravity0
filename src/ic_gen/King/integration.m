P.k = 1;
P.j = 1;
P.V0 = -1;
P.G = 1;

rspan = linspace(0.01, 1.2, 100);

q0 = [P.V0; rspan(1)*4*pi*P.G*rho(P.V0,P)/3];

opts = [];

[t,q] = ode45(@(t,q) f(t,q,P),rspan,q0,opts);

plot(rspan, q(:,1))
