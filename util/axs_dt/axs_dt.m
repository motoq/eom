function y = axs_dt(x)
%
% Input:
%   x  angular velocity, deg/minute
%
% Return:
%   Time step, minutes
%

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

x_dpm = (flip(theta_dot_deg_min))';
y_dtm = (flip(dt_minutes))';
A = (1./x_dpm);
k = ((A'*A)^-1)*A'*y_dtm;
y_shift = y_dtm(3) - k/x_dpm(3);

y = k./x + y_shift;

