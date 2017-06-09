function M=Tx(A)
    a1=1;       b1=0;       c1=0;
    a2=0;       b2=cos(A);  c2=sin(A);
    a3=0;       b3=-sin(A); c3=cos(A);
    M=[a1,b1,c1;a2,b2,c2;a3,b3,c3];
end
function M=Tz(A)
    a1=cos(A);  b1=sin(A);  c1=0;
    a2=-sin(A); b2=cos(A);  c2=0;
    a3=0;       b3=0;       c3=1;
    M=[a1,b1,c1;a2,b2,c2;a3,b3,c3];
end
function M=Ty(A)
    a1=cos(A);  b1=0; c1=-sin(A);
    a2=0;       b2=1; c2=0;
    a3=sin(A);  b3=0; c3=cos(A);
    M=[a1,b1,c1;a2,b2,c2;a3,b3,c3];
end