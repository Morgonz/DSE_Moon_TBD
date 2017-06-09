N=1000;
N0=1;
Y=zeros(N0,N);
X=zeros(1,N);

% for i=(1:N)
%     %    plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
%     Y(i)=plane_change(360/N*i,90,30,10000,700,0,90).run().incr;
%     X(i)=360/N*i;
% end

for i=(1:N)
    %    plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
    Y(i)=plane_change(90/N*i,90,30,10000,1736+700,0,90).run().incr;
    X(i)=90/N*i;
end

for i=(1:N)
    %    plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
    Y2(i)=plane_change(0,90/N*i,30,10000,1736+700,0,90).run().incr;
    X2(i)=90/N*i;
end

plot(X,Y); hold on; plot(X2,Y2);