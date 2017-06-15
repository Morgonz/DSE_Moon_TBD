function C=MPPVfun(ra)

    nplanes=3;           %number of planes to insert into
    %rleo=6671;          %starting LEO radius
    %ra=1737+1629;       %apoapse after initial manoeuvre
    rp=1737+1629;        %final moon orbit radius
    reol=1737+100;
    vinf=0.830;
    MuM=4.905*10^12;
    Mdry_sat=13.75;

    Isp=330;
    Isp_sat=220;
    Mbase=504;
    m=8;
    g=9.81;


    %constellation info
    INCC=50.2;
    AAND=360/6;
    A=zeros(nplanes,2);

    for i=(1:nplanes)
        A(i,1)=AAND*(i-1)+314.3418;
        A(i,2)=INCC;
    end

    w_p=70;

    dv2=zeros(nplanes-1,1);
    
    %dvs
    vp1=((vinf*10^3)^2+2*MuM/rp*10^(-3))^0.5;
    vp2=(MuM*(2/rp*10^(-3)-2/(ra+rp)*10^(-3)))^0.5;
    vc=(MuM/rp*10^(-3))^0.5;
    dv1=vp1-vp2;
    dv3=vp2-vc;
    va_eol=(MuM*(2/ra*10^(-3)-2/(ra+reol)*10^(-3)))^0.5;
    va=(MuM*(2/ra*10^(-3)-2/(ra+rp)*10^(-3)))^0.5;
    dv_eol=va-va_eol;
    
    %half plane change
    dv_hpc=489;
    Mrat_hpc=exp(dv_hpc/Isp/g);
    
    %sat mass
    Mrat_sat=exp((dv3)/Isp_sat/g);
    Msat=Mrat_sat*Mdry_sat;

    %PLANE CHANGES
    %plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
    F=zeros(nplanes-1,3);
    M=zeros(nplanes-1,3);
    for i=(1:(nplanes-1))
        plane=plane_change(A(i,1),A(i,2),w_p,ra,rp,A(i+1,1),A(i+1,2)).run();
        w_p=plane.aopF*360/2/pi;
        dv2(i)=plane.dv;
        F(i,:)=plane.F;
        M(i,:)=plane.m_pos;
        plane.cplot();
    end

    %vehicle mass estimation
    Mrat2=exp(dv2/Isp/g);
    Mrat_eol=exp(dv_eol/Isp/g);
    Mrat_ins=exp(dv1/Isp/g);
    Mt=Mbase*Mrat_eol;
    for i=(1:nplanes-1)
        Mt=(Mt+m*Msat)*Mrat2(i);
    end
    Mt=Mt*Mrat_hpc;
    Mt=Mt*Mrat_ins;
    C=F;
end
