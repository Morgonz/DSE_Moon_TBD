classdef getkepler
    properties
        mu=4.902801e12;
        rx;ry;rz;
        vx;vy;vz;
        
        hx;hy;hz;
        SMA;ECC;RAAN;INC;S;
    end
    methods
        
        function obj=getkepler(rx,ry,rz,vx,vy,vz)
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
            r=sqrt(obj.rx.^2+obj.ry.^2+obj.rz.^2);
            e1=(((obj.vy.*obj.hz)-(obj.vz.*obj.hy))/obj.mu)-(obj.rx./r);
            e2=(((obj.vz.*obj.hx)-(obj.vx.*obj.hz))/obj.mu)-(obj.ry./r);
            e3=(((obj.vx.*obj.hy)-(obj.vy.*obj.hx))/obj.mu)-(obj.rz./r);
            obj.ECC=sqrt(e1.^2+e2.^2+e3.^2);
            
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
        end
        
    end
end