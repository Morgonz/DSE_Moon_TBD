function an=tanomaly(an0,dt,sma,e,tol)
mu=3.9860044e14;

n=(mu/sma^3)^.5;

Ein=2*atan(((1-e)/(1+e))^.5*tan(an0/2));
tin=1/n*(Ein-e*sin(Ein));
t=tin+dt;

M=n*t-floor(n*t/2/pi)*2*pi;

if M>pi
    E0=M+e/2;
else
    E0=M-e/2;
end

%run=true;
% while run
%     DE0=-(E0-e*sin(E0)-n*t)/(1-e*cos(E0));
%     E1=E0+DE0;
%     if DE0<tol
%         run=false;
%     end
% end
run=true;
while run
    E1=e*sin(E0)+M;
    DE0=E1-E0;
    E0=E1;
    if DE0<tol
        run=false;
    end
end

an=2*atan(((1+e)/(1-e))^.5*tan(E1/2));

if M>pi
    an=2*pi+an;
end

end

%tanomaly(0.1,5,(6371+300)*10^3,0.045,0.001)