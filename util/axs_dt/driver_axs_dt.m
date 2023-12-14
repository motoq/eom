% Script to visualize and fit an appropriate time increment to use when
% searching for an orbit to ground point access window given the rate of
% change in true anomaly
%
% Kurt Motekew  2023/12/12
%

close all;
clear;

  % 8 seconds at low LEO, 2 minutes at Geo
dt_minutes = [8 16 120 240]/60;
p_minutes = [90 150 720 1440];

re_km = 6378.1370;
mu_km3_s2 = 398600.4415;
p_sec = 60*p_minutes;
a3_km = mu_km3_s2*(p_sec.*p_sec)/(4*pi*pi);
a_km = nthroot(a3_km, 3);
a_re = a_km/re_km;

v_cir_km_s = sqrt(mu_km3_s2./a_km);
theta_dot_sec = v_cir_km_s./a_km;
theta_dot_deg_min = 60*180*theta_dot_sec/pi;

  % Get a feel for the data
figure;  hold on;
plot(p_minutes, dt_minutes, 'k-o',...
     p_minutes, a_re, 'b-o',...
     p_minutes, theta_dot_deg_min, 'r-o');
xlabel('Period (minutes)');
ylabel('dt (minutes)');
title('dt vs. Period');
legend('dt','a (ER)', 'deg/min');

  % Fit time increment vs. angular rate
x_dpm = (flip(theta_dot_deg_min))';
y_dtm = (flip(dt_minutes))';
A = (1./x_dpm);
k = ((A'*A)^-1)*A'*y_dtm;
y_shift = y_dtm(3) - k/x_dpm(3);
x = x_dpm(1):.1:x_dpm(end);
y = k./x + y_shift;
  % and plot
figure;  hold on;
plot(x_dpm, y_dtm, 'k-o');
plot(x, y, 'r-');
xlabel('dnu/dt (deg/minute)');
ylabel('deltaT (minutes)');
title('Access Search Time Increment vs. True Anomaly Rate');
legend('Data', 'Fit');

fprintf('\nminutes = %1.5f/(deg/minute) + %1.4f', k, y_shift);

  % Plot with time increment in seconds
figure;  hold on;
plot(x, 60*y, 'r-');
xlabel('dnu/dt (deg/minute)');
ylabel('deltaT (seconds)');
title('Access Search Time Increment vs. True Anomaly Rate');

fprintf('\n\n');
