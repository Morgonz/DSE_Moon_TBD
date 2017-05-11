function A = AfromSymbolic(rx,ry,rz,dotrx,dotry,dotrz, rmu)
A(1,1) = 0;
A(1,2) = 0;
A(1,3) = 0;
A(1,4) = 1;
A(1,5) = 0;
A(1,6) = 0;


A(2,1) = 0;
A(2,2) = 0;
A(2,3) = 0;
A(2,4) = 0;
A(2,5) = 1;
A(2,6) = 0;


A(3,1) = 0;
A(3,2) = 0;
A(3,3) = 0;
A(3,4) = 0;
A(3,5) = 0;
A(3,6) = 1;


A(4,1) = (rmu - 1)/((rmu + rx)^2 + ry^2 + rz^2)^(3/2) - rmu/((rmu + rx - 1)^2 + ry^2 + rz^2)^(3/2) + (3*rmu*(2*rmu + 2*rx - 2)*(rmu + rx - 1))/(2*((rmu + rx - 1)^2 + ry^2 + rz^2)^(5/2)) - (3*(2*rmu + 2*rx)*(rmu + rx)*(rmu - 1))/(2*((rmu + rx)^2 + ry^2 + rz^2)^(5/2)) + 1;
A(4,2) = (3*rmu*ry*(rmu + rx - 1))/((rmu + rx - 1)^2 + ry^2 + rz^2)^(5/2) - (3*ry*(rmu + rx)*(rmu - 1))/((rmu + rx)^2 + ry^2 + rz^2)^(5/2);
A(4,3) = (3*rmu*rz*(rmu + rx - 1))/((rmu + rx - 1)^2 + ry^2 + rz^2)^(5/2) - (3*rz*(rmu + rx)*(rmu - 1))/((rmu + rx)^2 + ry^2 + rz^2)^(5/2);
A(4,4) = 0;
A(4,5) = 2;
A(4,6) = 0;


A(5,1) = (3*rmu*ry*(2*rmu + 2*rx - 2))/(2*((rmu + rx - 1)^2 + ry^2 + rz^2)^(5/2)) - (3*ry*(2*rmu + 2*rx)*(rmu - 1))/(2*((rmu + rx)^2 + ry^2 + rz^2)^(5/2));
A(5,2) = (rmu - 1)/((rmu + rx)^2 + ry^2 + rz^2)^(3/2) - rmu/((rmu + rx - 1)^2 + ry^2 + rz^2)^(3/2) - (3*ry^2*(rmu - 1))/((rmu + rx)^2 + ry^2 + rz^2)^(5/2) + (3*rmu*ry^2)/((rmu + rx - 1)^2 + ry^2 + rz^2)^(5/2) + 1;
A(5,3) = (3*rmu*ry*rz)/((rmu + rx - 1)^2 + ry^2 + rz^2)^(5/2) - (3*ry*rz*(rmu - 1))/((rmu + rx)^2 + ry^2 + rz^2)^(5/2);
A(5,4) = -2;
A(5,5) = 0;
A(5,6) = 0;


A(6,1) = (3*rmu*rz*(2*rmu + 2*rx - 2))/(2*((rmu + rx - 1)^2 + ry^2 + rz^2)^(5/2)) - (3*rz*(2*rmu + 2*rx)*(rmu - 1))/(2*((rmu + rx)^2 + ry^2 + rz^2)^(5/2));
A(6,2) = (3*rmu*ry*rz)/((rmu + rx - 1)^2 + ry^2 + rz^2)^(5/2) - (3*ry*rz*(rmu - 1))/((rmu + rx)^2 + ry^2 + rz^2)^(5/2);
A(6,3) = (rmu - 1)/((rmu + rx)^2 + ry^2 + rz^2)^(3/2) - rmu/((rmu + rx - 1)^2 + ry^2 + rz^2)^(3/2) - (3*rz^2*(rmu - 1))/((rmu + rx)^2 + ry^2 + rz^2)^(5/2) + (3*rmu*rz^2)/((rmu + rx - 1)^2 + ry^2 + rz^2)^(5/2);
A(6,4) = 0;
A(6,5) = 0;
A(6,6) = 0;

