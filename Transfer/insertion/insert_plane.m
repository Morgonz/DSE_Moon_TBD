%in
r_leo=6671;     %initial orbit about the Earth
r_pea=1880;     %PeA about the moon where second manoeuvre takes place
r_apa=9693;     %ApA about the moon after second manoevre is made
AAN_F=0;        %intended final Argument of Ascending Node
INC_F=45;       %intended final inclination

%guess
phi_int=20;     %intercept angle
W_int=30;       %intercept longitude


transfer=direct0(r_leo,23.4,r_pea).run();

insertion=insertion_coords(W_int,phi_int,r_pea,transfer.vinf).run();

plane=plane_change(insertion.W_an,insertion.inc,insertion.w_pea,r_apa,r_pea,AAN_F,INC_F).run();

plane.dv