classdef direct0
    properties
        r0;     %initial orbit
        inc0;   %initial inclination (ecliptic)
        rf;     %final orbit
        %----out
        dv1;
        dv2;
        dvt;
        inc_inf_r;
        vinf;
        debug;
    end
    methods
        function obj=direct0(a,b,c) %   initialise (r0,inc0,rf)
            obj.r0=a;
            obj.inc0=b;
            obj.rf=c;
        end
        function obj=run(obj)
            % CR3B EMTransfer arriving at node

            % define prking orbit
            %   r0 [km] 
            %   inc [deg] (ecliptic)
            % define final orbit
            %   rf [km]
            obj.r0=obj.r0*10^3;
            obj.rf=obj.rf*10^3;
            obj.inc0=obj.inc0*2*pi/360;

            %orbit parameters
            MuE=3.986*10^14;
            MuM=4.905*10^12;
            incm=5.145/360*2*pi;
            inc=obj.inc0-incm;
            rm=384399*10^3;
            vm=(MuE/rm)^(1/2);
            aT=(obj.r0+rm)/2;
            %velocities
            v0=(MuE/obj.r0)^(1/2);
            vT1=(MuE*(2/obj.r0-1/aT))^(1/2);
            vT2=(MuE*(2/rm-1/aT))^(1/2);
            obj.vinf=((vm-vT2*cos(inc))^2+(vT2*sin(inc))^2)^(1/2);
            obj.inc_inf_r=asin(vT2*sin(inc)/(vm-vT2*cos(inc)));
            vf1=(obj.vinf^2+2*MuM/obj.rf)^(1/2);
            vf2=(MuM/obj.rf)^(1/2);
            %DVs
            obj.dv1=vT1-v0;
            obj.dv2=vf1-vf2;
            obj.dvt=obj.dv1+obj.dv2;
            obj.debug=0;
            obj.inc_inf_r=obj.inc_inf_r/2/pi*360;
        end
    end
end

    %a=direct0(6671,28.58,1736+700).run();