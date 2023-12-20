function orbital_w(a, e)
% ORBITAL_W is used to study a quick way of determining what time
% increment should be used to search for an orbit to ground point
% access window given only the radial distance of the satellite.
%
% Given an orbit, perigee and apogee distances, and angular rates, a
% simple linear interpoloation function of radial distance is created
% that returns an appropriate search interval time increment.  This is a
% cheap yet accurate way of computing the time increment given the
% satellite position vector needs to be computed at each evaluation time
% anyway.  In contrast, actually computing the angular velocity w.r.t.
% the center of the earth requires pulling both position and velocity,
% along with performing a number of other computations when the goal is
% to minimize them.
%
% Inputs:
%   a  Semimajor axis, ER
%   e  Eccentricity
%
% Kurt Motekew  2023/12/12
%

    % Constants and conversions - work in natural orbital units
    % but output in human readable ones
  gm = 1.0;
  sec_per_tu = 806.811;
  tu_per_sec = 1/sec_per_tu;
  tu_per_min = 60*tu_per_sec;
  rad_per_deg = pi/180;

    % Additional orbital parameters
  b = a*sqrt(1 - e*e);
  phi = 0;

    % Interpolation parameters
  rp = a*(1 - e);
  ra = a*(1 + e);
  vp = sqrt(gm*(2/rp - 1/a));
  va = sqrt(gm*(2/ra - 1/a));
  theta_dot_p = vp/rp;
  theta_dot_a = va/ra;

    % ...vs. true anomaly
  p = a*(1 - e*e);
  nu_s = pi*(0:2:360)/180;
  r_s = p./(1 + e*cos(nu_s));
  v_s = sqrt(gm*(2./r_s - 1/a));
  den_s = sqrt(1 + 2*e.*cos(nu_s) + e*e);
    % flight path angle
  fpa_s = atan2(e.*sin(nu_s)./den_s, (1 + e.*cos(nu_s))./den_s);

    % This is what we would want to typically compute as the
    % independent variable leading to the time increment
  theta_dot_s = (v_s.*cos(fpa_s))./r_s;

    % dt vs. theta_dot equation with upper and lower bounds
  k = 0.99547*rad_per_deg;
  c = -0.1481*tu_per_min;
  lb = 8.0*tu_per_sec;
  ub = 2.0*tu_per_min;

    % Full algorithm based dt
  dt = min(max(k./theta_dot_s + c, lb), ub);

    % Linear interpolation based dt
  dt_a = min(max(k./theta_dot_a + c, lb), ub);
  dt_p = min(max(k./theta_dot_p + c, lb), ub);
  sf = (r_s - rp)/(ra - rp);
  dt_approx = dt_p + sf*(dt_a - dt_p);

    % Get semilatus rectum theta_dot
  fpa_slr = atan2(e/sqrt(1 + e*e), 1/sqrt(1 + e*e));
  r_slr = p;
  v_slr = sqrt(gm*(2/r_slr - 1/a));
  v_slr_t = v_slr*cos(fpa_slr);
  theta_dot_slr = v_slr_t/r_slr;
  dt_slr = min(max(k/theta_dot_slr + c, lb), ub);
  Ap = [1 rp 10^rp; 1 r_slr 10^r_slr ; 1 ra 10^ra];
  y = [dt_p ; dt_slr ; dt_a];
  W = [1 0 0 ; 0 1 0 ; 0 0 1];
  phat = (Ap'*W*Ap)^-1*(Ap'*y);
  dt_approx_2 = min(max(phat(1) + phat(2).*r_s + phat(3).*(10.^r_s), lb), ub);

    % Plot orbit shape
  figure;  hold on;
  scatter(0, 0, 'b');
  plot(r_s.*cos(nu_s), r_s.*sin(nu_s), 'k-');
  plot(cos(nu_s), sin(nu_s), 'b-');
  xlabel('x');
  ylabel('y');
  title('Orbit');
  axis equal;

  figure;  hold on;
  plot(180*nu_s/pi, 180*fpa_s/pi);
  xlabel('(deg)');
  xlabel('(deg)');
  title('Flight Path Angle vs. True Anomaly');
  xlim([0 360]);

  figure;  hold on;
  plot(180*nu_s/pi, tu_per_sec*180*theta_dot_s/pi, 'b-');
  scatter(0, tu_per_sec*180*theta_dot_p/pi, 'b');
  scatter(180, tu_per_sec*180*theta_dot_a/pi, 'b');
  scatter(360, tu_per_sec*180*theta_dot_p/pi, 'b');
  scatter(90, tu_per_sec*180*theta_dot_slr/pi, 'm');
  scatter(270, tu_per_sec*180*theta_dot_slr/pi, 'm');
  xlabel('(deg)');
  ylabel('(deg/sec)');
  title('True Anomaly Dot vs. True Anomaly');
  xlim([0 360]);

    % Full eqn vs. interpolated
  figure;  hold on;
  plot(180*nu_s/pi, sec_per_tu*dt, 'b-');
  plot(180*nu_s/pi, sec_per_tu*dt_approx, 'r-');
  plot(180*nu_s/pi, sec_per_tu*dt_approx_2, 'm-');
  xlabel('(deg)');
  ylabel('(sec)');
  title('Search interval Increment vs. True Anomaly');
  xlim([0 360]);
  legend('Full Algorithm', 'Interpolation', '10^x Interpolation');

