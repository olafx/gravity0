function result = rho(V, P)

k = P.k;
j = P.j;
V0 = P.V0;

result = pi^1.5*k/(j^3)*exp(2*j^2*(V0-V))*erf(j*sqrt(-2*V))-2*pi*k*sqrt(-2*V)*exp(2*j^2*V0)*(1/(j^2)-4*V/3);

end
