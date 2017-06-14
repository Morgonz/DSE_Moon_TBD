classdef getkepler_e
    properties
        mu=4.0350323550225975e14;
        rx;ry;rz;
        vx;vy;vz;
        
        hx;hy;hz;
        SMA;ECC;RAAN;INC;AOP;S;R;nor;ta;evecx;evecy;evecz;
    end
    methods
        
        function obj=getkepler_e(rx,ry,rz,vx,vy,vz)
            siz=size(rx);siz=siz(1);
            obj.S=siz;
            obj.rx=rx; obj.ry=ry; obj.rz=rz;
            obj.vx=vx; obj.vy=vy; obj.vz=vz;
            
            
            %ANGULAR MOMENTUM
            obj.hx=(obj.ry.*obj.vz)-(obj.rz.*obj.vy);
            obj.hy=(obj.rz.*obj.vx)-(obj.rx.*obj.vz);
            obj.hz=(obj.rx.*obj.vy)-(obj.ry.*obj.vx);
            
            %SEMI MAJOR AXIS
            obj.SMA=zeros(obj.S,1);
            for i=(1:obj.S)
                af1=1/(obj.rx(i)^2+obj.ry(i)^2+obj.rz(i)^2)^.5*2;
                af2=(obj.vx(i)^2+obj.vy(i)^2+obj.vz(i)^2)/obj.mu;
                obj.SMA(i)=1/(af1-af2);
            end
            
            %ECCENTRICITY
            obj.ECC=zeros(obj.S,1);
            obj.R=sqrt(obj.rx.^2+obj.ry.^2+obj.rz.^2);
            ex=(((obj.vy.*obj.hz)-(obj.vz.*obj.hy))/obj.mu)-(obj.rx./obj.R);
            ey=(((obj.vz.*obj.hx)-(obj.vx.*obj.hz))/obj.mu)-(obj.ry./obj.R);
            ez=(((obj.vx.*obj.hy)-(obj.vy.*obj.hx))/obj.mu)-(obj.rz./obj.R);
            obj.ECC=sqrt(ex.^2+ey.^2+ez.^2);
            
            obj.evecx=ex;
            obj.evecy=ey;
            obj.evecz=ez;
            
            %INCLINATION
            obj.INC=(acos(obj.hz./(sqrt(obj.hx.^2+obj.hy.^2+obj.hz.^2))));
            
            %ASCENDING NODE
            nx=zeros(obj.S,1);
            ny=zeros(obj.S,1);
            for i=(1:obj.S)
                nx(i)=-obj.hy(i);
                ny(i)=obj.hx(i);
            end
            
            %RAAN
            obj.RAAN=zeros(obj.S,1);
            for i=(1:obj.S)
                obj.RAAN(i)=acos(nx(i)/(nx(i)^2+ny(i)^2)^.5);
                if ny(i)<0
                    obj.RAAN(i)=2*pi-obj.RAAN(i);
                end
            end
            
            %ARGUMENT OF PERIAPSE
            obj.AOP=zeros(obj.S,1);
            for i=(1:obj.S)
                a=nx(i)*ex(i);
                b=ny(i)*ey(i);
                d=(ex(i)^2+ey(i)^2+ez(i)^2)*(nx(i)^2+ny(i)^2);
                obj.AOP(i)=acos((a+b)/d^.5);
                if ez(i)<0
                    obj.AOP(i)=2*pi-obj.AOP(i);
                end
            end
            
            %NORMAL VECTOR
            obj.nor=zeros(obj.S,3);
            for i=(1:obj.S)
                obj.nor(i,:)=cross([obj.rx(i);obj.ry(i);obj.rz(i)],[obj.vx(i);obj.vy(i);obj.vz(i)]);
                obj.nor(i,:)=obj.nor(i,:)/(obj.nor(i,1)^2+obj.nor(i,2)^2+obj.nor(i,3)^2)^.5;
            end
            
            %TRUE ANOMALY
            obj.ta=zeros(obj.S,1);
            for i=(1:obj.S)
                obj.ta(i)=acos((ex(i)*rx(i)+ey(i)*ry(i)+ez(i)*rz(i))/(obj.ECC(i)*obj.R(i)));
                if vz<0
                    obj.ta(i)=2*pi-obj.ta(i);
                end
            end
            
        end
        
    end
end