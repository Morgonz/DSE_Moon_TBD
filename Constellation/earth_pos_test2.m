T=60*60*24*27;
n=100;
Re=zeros(n,3);
tl=zeros(n,1);
re=zeros(n,1);
aban=zeros(n,1);
tan=zeros(n,1);
for i=(1:n)
    t=i/n*T;
    tl(i)=t;
    
    RAANe=1.4180;
    INCe=0.1177;
    theta_in=1.5835;
    an0=0.5988;
    e_m=0.0450;
    sma_m=3.8140e+08;

    an=tanomaly(an0,t,sma_m,e_m,.0001);
    theta=an-an0+theta_in;
    

    T3=rotationmatrices.Tz(RAANe);
    T2=rotationmatrices.Tx(INCe);
    T1=rotationmatrices.Tz(theta);

    Re(i,:)=T3*T2*T1*[1;0;0]*sma_m*(1-e_m^2)/(1+e_m*cos(theta));
    re(i)=(Re(i,1)^2+Re(i,2)^2+Re(i,3)^2)^.5;
    
    tan(i)=theta;
    aban(i)=acos(Re(i,3)/re(i));
end

%%

scatter3(Re(:,1),Re(:,2),Re(:,3))
axis equal;
axis vis3d;

%%
scatter(tl,Re(:,1))

%%
scatter(tl,re(:,1))

%%
scatter(tl,aban)

%%
scatter(tl,tan)