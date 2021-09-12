function result = f(r,q,P)

k = P.k;
j = P.j;
V0 = P.V0;
G = P.G;

result = [q(2)
          4*pi*G*rho(q(1),P)-2*q(2)/r];
end
