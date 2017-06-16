theta_deg = 0
theta_rad = theta_deg*(pi/180)
r = 1
r_dot = 0
theta_dot = 1


T = [ cos(theta_rad) -sin(theta_rad);
      sin(theta_rad) cos(theta_rad)  ]
  
%V_polar = [r_dot  ; r*theta_dot ]

%V_cartesian = T*V_polar
V_cartesian = [1; 0]

V_polar = inv(T)*V_cartesian
