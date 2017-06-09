classdef insertion_coords
    properties
        %----IN
        W_int;      %longitude where s/c intercepts Moon
        phi_int;    %intercept angle (in Vinf axis)
        r_pea;      %periapse radius (where manoeuvre is made)
        v_inf;      %assume vinf for hyperbolic trajectory
        %----OUT
        inc;        %inclination
        W_an;       %longitude of ascending node
        w_pea;      %argument of periapse
        D;          %distance between inf v vector and node axis
        %----constants
        MuM=4.905*10^12;
        
    end
    methods
        function obj=insertion_coords(W_int,phi_int,r_pea,v_inf)
            obj.W_int=W_int*2*pi/360;
            obj.phi_int=phi_int*2*pi/360;
            obj.r_pea=r_pea*10^3;
            obj.v_inf=v_inf*10^3;
        end
        function obj=run(obj)
            %check exceptions
            if abs(obj.phi_int)>pi/2        %-90<phi_int<90
                if obj.phi_int>0
                    obj.phi_int=pi/2;
                else
                    obj.phi_int=-pi/2;
                end
            end
            
            a=-obj.MuM/obj.v_inf^2;
            e=1-obj.r_pea/a;
            obj.D=-a*(e^2-1)^(1/2);
            theta_inf=acos(-1/e);
            
            if obj.phi_int<0
                obj.W_an=obj.W_int+pi/2;
                obj.w_pea=theta_inf-pi;
            else
                obj.W_an=obj.W_int;
                obj.w_pea=theta_inf;
            end
            obj.inc=abs(obj.phi_int);
            %change units
            obj.W_an=obj.W_an/2/pi*360;
            obj.w_pea=obj.w_pea/2/pi*360;
            obj.inc=obj.inc/2/pi*360;
        end
    end
end
%name=insertion_coords(10,30,1900,0.8438188).run()