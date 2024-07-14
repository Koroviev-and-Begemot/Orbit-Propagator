function [a,e,i,RAAN,nu,omega] = RandV2Elements(R,V,mu)

r = norm(R)                     ;  v = norm(V)                  ;

% Angular Momentum
H = cross(R,V);
h = norm(H);

%Semi Latus Rectum
p = h^2/r;

% SMA and eccentricity
a = mu/(2*mu/r - v^2);
E = 1/mu*(cross(V,H)- mu/r*R);
e = norm(E);

%inclination
i = acosd(dot(H,[0 0 1])/(h));

%Node line
N = cross([0 0 1],H);
n = norm(N);

%RAAN
if N(2) > 0
    RAAN = acosd(N(1)/norm(N));
else
    RAAN = 360 - acosd(N(1)/norm(N));
end

%True argument of latitude
if R(3) > 0
    nu = acosd(dot(N,R)/(n*r));
else
    nu = 360 - acosd(dot(N,R)/(n*r));
end

%Argument Perigee
if E(3) < 0
    omega = 360 - acosd(dot(N,E)/(n*e));
else
    omega = acosd(dot(N,E)/(n*e));
end

end

