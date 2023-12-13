function orbital_w(a, e)

  gm = 1.0;
  sec_per_tu = 806.811;
  tu_per_sec = 1/sec_per_tu;
  tu_per_min = 60*tu_per_sec;
  rad_per_deg = pi/180;

%  a = 2;
%  e = .8;
%  e = .7;
  b = a*sqrt(1 - e*e);
  phi = 0;

  rp = a*(1 - e);
  ra = a*(1 + e);
  vp = sqrt(gm*(2/rp - 1/a));
  va = sqrt(gm*(2/ra - 1/a));
  theta_dot_p = vp/rp;
  theta_dot_a = va/ra;
  p = a*(1 - e*e);

  nu_s = pi*(0:2:360)/180;
  r_s = p./(1 + e*cos(nu_s));
  v_s = sqrt(gm*(2./r_s - 1/a));
  den_s = sqrt(1 + 2*e.*cos(nu_s) + e*e);
  fpa_s = atan2(e.*sin(nu_s)./den_s, (1 + e.*cos(nu_s))./den_s);

  theta_dot_s = (v_s.*cos(fpa_s))./r_s;


  k = 0.99547*rad_per_deg;
  c = -0.1481*tu_per_min;
  lb = 8.0*tu_per_sec;
  ub = 2.0*tu_per_min;

  dt = min(max(k./theta_dot_s + c, lb), ub);

  dt_a = min(max(k./theta_dot_a + c, lb), ub);
  dt_p = min(max(k./theta_dot_p + c, lb), ub);
  sf = (r_s - rp)/(ra - rp);
  dt_approx = dt_p + sf*(dt_a - dt_p);


  



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
  xlabel('(deg)');
  ylabel('(deg/sec)');
  title('True Anomaly Dot vs. True Anomaly');
  xlim([0 360]);

  figure;  hold on;
  %plot(tu_per_min*180*theta_dot_s/pi, sec_per_tu*dt, 'b-');
  %plot(tu_per_min*180*theta_dot_s/pi, sec_per_tu*dt_approx, 'r-');
  plot(180*nu_s/pi, sec_per_tu*dt, 'b-');
  plot(180*nu_s/pi, sec_per_tu*dt_approx, 'r-');
  xlabel('(deg)');
  ylabel('(sec)');
  title('Search interval Increment vs. True Anomaly');

