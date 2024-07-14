clear; close all
%% Constants and initialization

mu_earth = 398600  ;
rad_earth = 6371   ;
mu_moon   = 4.904E3;
rad_moon  = 1737   ;

delta_t = 10;
% sim_time = 0.1*2.5E6;
sim_time = 31*24*60^2;
t = 0;
k = ceil(sim_time/delta_t)+1;

moon_state         = zeros(6,k);
moon_state_dot     = zeros(6,k);
moon_state_bar     = zeros(6,k);
moon_state_dot_bar = zeros(6,k);

sat_state         = zeros(6,k);
sat_state_dot     = zeros(6,k);
sat_state_bar     = zeros(6,k);
sat_state_dot_bar = zeros(6,k);
a_b               = zeros(3,k);

m_sat = 1000;
T_sat = 1;

t_b = [4*60^2 10*60^2];
tau_b = [500 5000];
dir = [1 2];
% 1: prograde 2: retrograde 3: normal 4: anti-normal

t_b = [t_b sim_time];
tau_b = [tau_b 0];
dir = [dir 1];

moon_state(:,1) = peri2ECI(384314.45,0.026,18.28,0,0,0,mu_earth);
moon_state_dot(:,1) = [moon_state(4:6,1); -mu_earth*moon_state(1:3,1)./norm(moon_state(1:3,1)).^3];

sat_state(:,1) = peri2ECI(35000,0.866,80,270,0,0,mu_moon) + moon_state(:,1);
%sat_state(:,1) = peri2ECI(7000,0,77,0,0,0,mu_earth);
sat_state_dot(:,1) = [sat_state(4:6,1); -mu_earth*sat_state(1:3,1)./norm(sat_state(1:3,1)).^3 - mu_moon*(sat_state(1:3,1) - moon_state(1:3,1))./norm(sat_state(1:3,1) - moon_state(1:3,1)).^3];

%% calculation loop

i = 2;
p = 1;

f = waitbar(0, 'Starting');

for j = 1:length(t_b)
    while t < t_b(j)

        moon_state_bar(:,i)     = moon_state(:,i-1) + delta_t*moon_state_dot(:,i-1);
        moon_state_dot_bar(:,i) = [moon_state_bar(4:6,i); -mu_earth*moon_state_bar(1:3,i)./norm(moon_state_bar(1:3,i)).^3];
        moon_state(:,i)         = moon_state(:,i-1) + delta_t/2*(moon_state_dot(:,i-1) + moon_state_dot_bar(:,i));
        moon_state_dot(:,i)     = [moon_state(4:6,i); -mu_earth*moon_state(1:3,i)./norm(moon_state(1:3,i)).^3];
    
        sat_state_bar(:,i)      = sat_state(:,i-1) + delta_t*sat_state_dot(:,i-1);
        sat_state_dot_bar(:,i)  = [sat_state_bar(4:6,i);
                                   - mu_earth*sat_state_bar(1:3,i)./norm(sat_state_bar(1:3,i)).^3 ...
                                   - mu_moon*(sat_state_bar(1:3,i) - moon_state_bar(1:3,i))./norm(sat_state_bar(1:3,i) ...
                                   - moon_state_bar(1:3,i)).^3];
        sat_state(:,i)          = sat_state(:,i-1) + delta_t/2*(sat_state_dot(:,i-1) + sat_state_dot_bar(:,i));
        sat_state_dot(:,i)      = [sat_state(4:6,i);
                                   - mu_earth*sat_state(1:3,i)./norm(sat_state(1:3,i)).^3 ...
                                   - mu_moon*(sat_state(1:3,i) ...
                                   - moon_state(1:3,i))./norm(sat_state(1:3,i) ...
                                   - moon_state(1:3,i)).^3];
        
        t = t + delta_t;
        i = i+1;


        if mod(i,1000) == 0
           out = ['Calculating ' num2str(round(i/k,2)*100) '%'];
           waitbar(i/k,f,out);
        end


    end
    
    while t < t_b(j) + tau_b(j)
        
        if tau_b(j) > 0
            a_b(:,i) = T_sat/m_sat/1000*sat_state(4:6,i-1)./(norm(sat_state(4:6,i-1)));
            switch dir(j)
                case 1
                    break
                case 2
                    a_b(:,i) = -a_b(:,i);
            end
        else
            a_b(:,i) = 0;
        end

        moon_state_bar(:,i)     = moon_state(:,i-1) + delta_t*moon_state_dot(:,i-1);
        moon_state_dot_bar(:,i) = [moon_state_bar(4:6,i); -mu_earth*moon_state_bar(1:3,i)./norm(moon_state_bar(1:3,i)).^3];
        moon_state(:,i)         = moon_state(:,i-1) + delta_t/2*(moon_state_dot(:,i-1) + moon_state_dot_bar(:,i));
        moon_state_dot(:,i)     = [moon_state(4:6,i); -mu_earth*moon_state(1:3,i)./norm(moon_state(1:3,i)).^3];
    
        sat_state_bar(:,i)      = sat_state(:,i-1) + delta_t*sat_state_dot(:,i-1);
        sat_state_dot_bar(:,i)  = [sat_state_bar(4:6,i);
                                   - mu_earth*sat_state_bar(1:3,i)./norm(sat_state_bar(1:3,i)).^3 ...
                                   - mu_moon*(sat_state_bar(1:3,i) ...
                                   - moon_state_bar(1:3,i))./norm(sat_state_bar(1:3,i) ...
                                   - moon_state_bar(1:3,i)).^3 ...
                                   + a_b(:,i)];
        sat_state(:,i)          = sat_state(:,i-1) + delta_t/2*(sat_state_dot(:,i-1) + sat_state_dot_bar(:,i));
        sat_state_dot(:,i)      = [sat_state(4:6,i); 
                                   - mu_earth*sat_state(1:3,i)./norm(sat_state(1:3,i)).^3 ...
                                   - mu_moon*(sat_state(1:3,i) ...
                                   - moon_state(1:3,i))./norm(sat_state(1:3,i) ...
                                   - moon_state(1:3,i)).^3 ...
                                   + a_b(:,i)];

        t = t + delta_t;
        i = i+1;


        if mod(i,1000) == 0
           out = ['Calculating ' num2str(round(i/k,2)*100) '%'];
           waitbar(i/k,f,out);
        end
    end
     disp(['burn ', num2str(j), ' occured at ', num2str(t), 's and lasted ', num2str(tau_b(j))])
end

close(f)
%% plots

figure(1)
hold on
grid on
set(gcf, 'Position',  [100, 100, 600, 600])
q = 4E5;
xlim([-1 1]*q)
ylim([-1 1]*q)
zlim([-1 1]*q)

drawsphere(0,0,0,rad_earth);
view(3)

 plot3(moon_state(1,:),moon_state(2,:),moon_state(3,:),'LineWidth',2,'Color','r')
 plot3(sat_state(1,:),sat_state(2,:),sat_state(3,:),'LineWidth',2,'Color','k')

%comet3(moon_state(1,:),moon_state(2,:),moon_state(3,:))
%comet3(sat_state(1,:),sat_state(2,:),sat_state(3,:))

figure(2)
hold on
grid on
set(gcf, 'Position',  [700, 100, 600, 600])
drawsphere(0,0,0,rad_moon);
q = 10E3;
xlim([-1 1]*q)
ylim([-1 1]*q)
zlim([-1 1]*q)

plot3(sat_state(1,:) - moon_state(1,:),sat_state(2,:) - moon_state(2,:),sat_state(3,:) - moon_state(3,:),'LineWidth',2,'Color','r')

for i = 1:length(t_b)-1
    xx = t_b(i)/delta_t:(t_b(i) + tau_b(i))/delta_t;
    plot3(sat_state(1,xx) - moon_state(1,xx),sat_state(2,xx) - moon_state(2,xx),sat_state(3,xx) - moon_state(3,xx),'LineWidth',2,'Color','b')
end

view(3)


function [out] = peri2ECI(a,e,i,omega,theta,OMEGA,mu)

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

out = [ri;vi];

end

function drawsphere(x,y,z,r)
    [a,b,c] = sphere;
    a = a*r + x;
    b = b*r + y;
    c = c*r + z;
    surf(a,b,c)
 end

