%define final orbit planes
OP=[0,0;0,0;0,0];   %OP=[aanF1,incF1;aanF2,incF2,... Plane1=>Plane2=>Plane3
%define maximum apa after initial insertion manoeuvre
apa_max=5000*10^3;
%constants
MuM=4.905*10^12;

%initialise variables
W_int=0;
inc_int=0;
size_op=size(OP);

for i=1:40
    for i=1:40
        in_plane=insertion_coords(10,30,1900,0.5).run();