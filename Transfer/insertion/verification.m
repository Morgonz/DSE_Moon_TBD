N=1000;
N0=6;
Y=zeros(N0,N);
X=zeros(1,N);

% for i=(1:N)
%     %    plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
%     Y(i)=plane_change(360/N*i,90,30,10000,700,0,90).run().incr;
%     X(i)=360/N*i;
% end

for i0=(1:N0)
    for i=(1:N)
        %    plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
        Y(i0,i)=plane_change(0,90/N0*i0,30,10000,1736+700,90/N*i,45).run().incr;
        X(i)=90/N*i;
    end
end



for i=(1:N0)
    plot(X,Y);
end
xlabel('\Delta \Omega_{an}  (degrees)')
ylabel('Relative inclination (degrees)')
title('Relative Inclination as a Function of \Omega_{an} (i_F=45)')
legend('i_I=15','i_I=30','i_I=45','i_I=60','i_I=75','i_I=90')

