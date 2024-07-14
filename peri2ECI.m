function [ri,vi] = peri2ECI(a,e,i,omega,theta,OMEGA,mu)

i = deg2rad(i);
OMEGA = deg2rad(OMEGA);
omega = deg2rad(omega);
theta = deg2rad(theta);

h = sqrt(a*mu*(1 - e^2));

rp = h^2/mu*1/(1 + e*cos(theta)) * [cos(theta); sin(theta); 0];
vp = mu/h * [-sin(theta); e + cos(theta); 0];

a1 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1]; a2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]; a3 = [cos(OMEGA) sin(OMEGA) 0; -sin(OMEGA) cos(OMEGA) 0; 0 0 1];

A = (a1 * a2 * a3)';

ri = A*rp;

vi = A*vp;

end

