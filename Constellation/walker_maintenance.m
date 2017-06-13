classdef walker_maintenance
    properties
        %moon parameters
        rm=1738100; mu=4.902801e12;
        %constellation parameters
        INC_c; n_planes; SMA_c; T;
        vc;
        %sim options
        options=odeset('RelTol',1e-6);
        %initial cartesian
        s0;                         %6x(n_planes) array
        %simulation results (cartesian)
        Rx;Ry;Rz;Vx;Vy;Vz;          %(n_planes)x1 cell
        %simulation resutls (keplerian)
        sma; ecc; inc; raan; aop;   %(n_planes)x1 cell
        %simulation time steps
        t;                          %(n_planes)x1 cell
        %design outputs
        dv;                         %(n_planes)x1 array
    end
    
    methods
        function obj=walker_maintenance(INC,n_planes,SMA,T)
            obj.INC_c=INC; obj.n_planes=n_planes; obj.SMA_c=SMA; obj.T=T;
            obj.vc=(obj.mu/SMA)^.5;
            
            %set initial conditions for each plane in simulation
            DRAAN=2*pi/obj.n_planes;
            obj.s0=zeros(6,obj.n_planes);
            for i=(1:obj.n_planes)
                s13=rotationmatrices.Tz(-DRAAN*i)*[0;cos(INC)*SMA;sin(INC)*SMA];
                s46=rotationmatrices.Tz(-DRAAN*i)*[obj.vc;0;0];
                obj.s0(:,i)=[s13(1);s13(2);s13(3);s46(1);s46(2);s46(3)];
            end
            
            %init output structure
            obj.t=cell(obj.n_planes,1);
            obj.Rx=cell(obj.n_planes,1);
            obj.Ry=cell(obj.n_planes,1);
            obj.Rz=cell(obj.n_planes,1);
            obj.Vx=cell(obj.n_planes,1);
            obj.Vy=cell(obj.n_planes,1);
            obj.Vz=cell(obj.n_planes,1);
            obj.sma=cell(obj.n_planes,1);
            obj.ecc=cell(obj.n_planes,1);
            obj.inc=cell(obj.n_planes,1);
            obj.raan=cell(obj.n_planes,1);
            obj.aop=cell(obj.n_planes,1);
            obj.dv=zeros(obj.n_planes,1);
            
        end
        
        function obj=prop(obj,plane)
            
            [ts,S] = ode45(@sat2BPdegH,[0,obj.T],obj.s0(:,plane),obj.options);
            rx=S(:,1); ry=S(:,2); rz=S(:,3);
            vx=S(:,4); vy=S(:,5); vz=S(:,6);
            
            obj.t{plane}=ts;
            
            obj.Rx{plane}=rx; obj.Ry{plane}=ry; obj.Rz{plane}=rz;
            obj.Vx{plane}=vx; obj.Vy{plane}=vy; obj.Vz{plane}=vz;
            
            kepler=getkepler(rx,ry,rz,vx,vy,vz);
            obj.sma{plane}=kepler.SMA;
            obj.ecc{plane}=kepler.ECC;
            obj.inc{plane}=kepler.INC;
            obj.raan{plane}=kepler.RAAN;
            obj.aop{plane}=kepler.AOP;
            
            obj=obj.maint_dv(plane);
            
            plane
        end
        
        function obj=prop_all(obj)
            for i=(1:obj.n_planes)
                obj=obj.prop(i);
            end
        end
        
        function obj=maint_dv(obj,plane)
            %final elements
            SMA1=obj.sma{plane}(end);
            ECC1=obj.ecc{plane}(end);
            INC1=obj.inc{plane}(end);
            RAAN1=obj.raan{plane}(end);
            AOP1=obj.raan{plane}(end);
            
            Apa1=SMA1*(1+ECC1);
            Pea1=SMA1*(1-ECC1);
            
            %target elements
            SMA0=obj.sma{plane}(1);
            ECC0=obj.ecc{plane}(1);
            INC0=obj.inc{plane}(1);
            RAAN0=obj.raan{plane}(1);
            AOP0=obj.raan{plane}(1);
            
            Apa0=SMA0*(1+ECC0);
            Pea0=SMA0*(1-ECC0);
            
            PCM=plane_change(RAAN1,INC1,AOP1,Apa1,Pea1,RAAN0,INC0).run();
            
            DVp=abs((obj.mu*(2/Pea1-2/(Apa1+Pea1)))^.5-(obj.mu*(2/Pea1-2/(Apa0+Pea1)))^.5);
            DVa=abs((obj.mu*(2/Apa0-2/(Apa0+Pea1)))^.5-(obj.mu*(2/Apa0-2/(Apa0+Pea0)))^.5);
            obj.dv(plane)=PCM.dv+DVp+DVa;
            
            
        end
        
    end
    
end

%a=walker_maintenance(50.2*pi/180,6,(1629+1737)*10^3,60*60*60*24);
    