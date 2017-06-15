classdef plane_change
    properties
        %----IN
        aanI; incI; apaI; peaI; aopI; aanF; incF; progI=true; progF=true;
        %----OUT
        rm; dv; theta1, theta2; incr; vpea; vtheta; aopF; AM; AMs; n_I=[1,0,0]; n_F=[0,1,0]; nF0; F; incr_deb; m_pos;
        %----Exceptions
        ex_same_plane; ex_opposite; ex_domain=0;
        %----other
        standard=true;
        %----const
        MuM=4.905*10^12;
    end
    methods (Static)
        %Rotation matrices
        function M=Tx(A)
            a1=1;       b1=0;       c1=0;
            a2=0;       b2=cos(A);  c2=sin(A);
            a3=0;       b3=-sin(A); c3=cos(A);
            M=[a1,b1,c1;a2,b2,c2;a3,b3,c3];
        end
        function M=Tz(A)
            a1=cos(A);  b1=sin(A);  c1=0;
            a2=-sin(A); b2=cos(A);  c2=0;
            a3=0;       b3=0;       c3=1;
            M=[a1,b1,c1;a2,b2,c2;a3,b3,c3];
        end
        function M=Ty(A)
            a1=cos(A);  b1=0; c1=-sin(A);
            a2=0;       b2=1; c2=0;
            a3=sin(A);  b3=0; c3=cos(A);
            M=[a1,b1,c1;a2,b2,c2;a3,b3,c3];
        end
    end
    methods
        function obj=plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
            if obj.standard
                obj.aanI=aanI/360*2*pi; obj.incI=incI/360*2*pi;
                obj.apaI=apaI*10^3; obj.aopI=aopI/360*2*pi;
                obj.peaI=peaI*10^3; obj.aanF=aanF/360*2*pi; obj.incF=incF/360*2*pi;
            end
        end
        function A=transform(obj,inc,aan,w,x,y,z)
            A=[x;y;z];
            A=obj.Tz(-aan)*obj.Tx(-inc)*obj.Tz(-w)*A;
%             aan
%             inc
        end
        function obj=run(obj)
            
            %check exceptions
            
            if obj.incI > pi          %0<inc<180 (inclination domain)
                obj.ex_domain=1;
            elseif obj.incI < 0
                obj.ex_domain=1;
            end
            if obj.incF > pi
                obj.ex_domain=1;
            elseif obj.incF < 0
                obj.ex_domain=1;
            end
            %check if planes are the same
            if obj.incI==obj.incF
                if obj.incI==0
                    obj.ex_same_plane=true;
                    if obj.progI && obj.progF
                        obj.ex_opposite=false;
                    else
                        obj.ex_opposite=true;
                    end
                elseif (obj.aanI==obj.aanF)
                    obj.ex_same_plane=true;
                    obj.ex_opposite=false;
                elseif (obj.aanI==obj.aanF+pi) || (obj.aanI==obj.aanF-pi)
                    obj.ex_same_plane=true;
                    obj.ex_opposite=true;   %check if orbits are opposite
                else
                    obj.ex_same_plane=false;
                    obj.ex_opposite=false;
                end
            else
                obj.ex_opposite=false;
                obj.ex_same_plane=false;
            end
            
            %intermediate orbit eccentricity
            e1=(obj.apaI-obj.peaI)/(obj.apaI+obj.peaI);
            
            %transformations
            D_OM=obj.aanI-obj.aanF;
            Tx_phi=obj.Tx(obj.incI);
            Tz_w=obj.Tz(obj.aopI);
            Tz_OM=obj.Tz(D_OM);
            obj.nF0=[0;-sin(obj.incF);cos(obj.incF)];
            nF=Tz_w*Tx_phi*Tz_OM*obj.nF0;
            
            %node angles
            obj.theta1=atan(abs(nF(1)/nF(2)));
            obj.theta2=pi-obj.theta1;
            
            %relative inclination
            obj.incr=acos(abs(nF(3)/(nF(1)^2+nF(2)^2+nF(3)^2)^(1/2)));
            obj.incr_deb=obj.incr;
            
            %ascending manoeuvre
            if ((nF(2)>0)&&(nF(1)<0))||((nF(2)<0)&&(nF(1)>0))
                obj.AM=false;
                obj.AMs=-1;
            else
                obj.AM=true;
                obj.AMs=1;
            end
            
            if nF(2)<0
                AMs2=1;
            else
                AMs2=-1;
            end
            
            %new argument of periapse
            MAI=obj.aopI-obj.AMs*obj.theta1;
            rf=obj.Tz(-D_OM)*obj.Tx(-obj.incI)*obj.Tz(-MAI)*[1;0;0];
            MAF=acos(rf(1));
                %MAF = acos(cos(AI)*cos(-D_OM)-sin(AI)*sin(D_OM)*(cos(-obj.incI)-sin(obj.incI)));
            obj.aopF=AMs2*MAF+obj.AMs*obj.theta1;
            
            %rinc displays the minimum angle
            if abs(obj.incr)>pi/2
                obj.incr=abs(pi-obj.incr);
            end
            
            if obj.ex_same_plane
                if obj.ex_opposite
                    vapa=(obj.MuM/obj.apaI)^(1/2);
                    obj.dv=2*vapa;      %dv to switch rotation
                    obj.rm=obj.apaI;
                else
                    obj.dv=0;
                    obj.rm=obj.apaI;
                end
            else
                obj.rm=((obj.peaI+obj.apaI)/2*(1-e1^2))/(1+e1*cos(obj.theta2));
                obj.vpea=(obj.MuM/obj.peaI)^(1/2);
                obj.vtheta=obj.vpea*obj.peaI/obj.rm;
                obj.dv=2*obj.vtheta*sin(obj.incr/2);
            end
            
            %Manoeuvre position
            xm=obj.rm*cos(obj.theta2); ym=obj.rm*obj.AMs*sin(obj.theta2);
            M=obj.transform(obj.incI,obj.aanI,obj.aopI,xm,ym,0);
            obj.m_pos=M;
            
            %Manoeuvre direction
            obj.n_F=obj.transform(0,obj.aanF,0,obj.nF0(1),obj.nF0(2),obj.nF0(3));
            nI0=[0;-sin(obj.incI);cos(obj.incI)];
            obj.n_I=obj.transform(0,obj.aanI,0,nI0(1),nI0(2),nI0(3));
            
            %R=cross(M/(M(1)^2+M(2)^2+M(3)^2)^.5,obj.n_I);
            R=cross(obj.n_I,obj.n_F);
            R=cross(R/(R(1)^2+R(2)^2+R(3)^2)^.5,obj.n_I);
            ratc=cos((pi-obj.incr)/2);
            rats=sin((pi-obj.incr)/2);
            obj.F=-ratc*R-rats*obj.n_I;
            
            
            
        end
        function obj=units2(obj)        %change units
            
            obj.incr=obj.incr*360/2/pi;
            obj.theta1=obj.theta1*360/2/pi;
            obj.theta2=obj.theta2*360/2/pi;
            obj.aopF=obj.aopF/2/pi*360;
        end
        function obj=cplot(obj)
            e1=(obj.apaI-obj.peaI)/(obj.apaI+obj.peaI);
            N=1000;
            X=zeros(1,N);
            Y=zeros(1,N);
            Z=zeros(1,N);
            X2=zeros(1,N);
            Y2=zeros(1,N);
            Z2=zeros(1,N);
            
            [P,Q,R]=sphere(40);
            P=P.*1736000;Q =Q.*1736000;R =R.*1736000;
            
            rm1=obj.rm*1.2;
            rm2=((obj.peaI+obj.apaI)/2*(1-e1^2))/(1+e1*cos(obj.theta1*pi/180))*1.2;
            Xn=[rm1*cos(obj.theta2),rm2*cos(obj.theta1)];
            Yn=[obj.AMs*rm1*sin(obj.theta2),-obj.AMs*rm2*sin(obj.theta1)];
            Zn=[0,0];
            for i=(1:2)
                M=obj.transform(obj.incI,obj.aanI,obj.aopI,Xn(i),Yn(i),Zn(i));
                Xn(i)=M(1);
                Yn(i)=M(2);
                Zn(i)=M(3);
            end
            Xn0=[rm1*cos(obj.theta2),rm2*cos(obj.theta1)];
            Yn0=[-obj.AMs*rm1*sin(obj.theta2),obj.AMs*rm2*sin(obj.theta1)];
            Zn0=[0,0];
            for i=(1:2)
                M=obj.transform(obj.incI,obj.aanI,obj.aopI,Xn0(i),Yn0(i),Zn0(i));
                Xn0(i)=M(1);
                Yn0(i)=M(2);
                Zn0(i)=M(3);
            end
            
            xm=obj.rm*cos(obj.theta2); ym=obj.rm*obj.AMs*sin(obj.theta2);
            M=obj.transform(obj.incI,obj.aanI,obj.aopI,xm,ym,0);
            xm=M(1); ym=M(2); zm=M(3);

            for i=(1:N)
                T=2*pi/N*i;
                r=((obj.peaI+obj.apaI)/2*(1-e1^2))/(1+e1*cos(T));
                x=r*cos(T);
                y=r*sin(T);
                z=0;
                M=obj.transform(obj.incI,obj.aanI,obj.aopI,x,y,z);
                X(i)=M(1);
                Y(i)=M(2);
                Z(i)=M(3);
                M=obj.transform(obj.incF,obj.aanF,obj.aopF,x,y,z);
                X2(i)=M(1);
                Y2(i)=M(2);
                Z2(i)=M(3);
            end
            
            nf=obj.n_F*3*10^6;
            ni=obj.n_I*3*10^6;
            fv=obj.F*2*10^6;
            
            %PLOT
            surf(P,Q,R);hold on;
            plot3(X,Y,Z);hold on;
            plot3(X2,Y2,Z2);hold on;
            
            plot3(Xn,Yn,Zn,'--','color','r');hold on;
            %plot3(Xn0,Yn0,Zn0,'--','color','g');hold on;
            
            %normal vectors
%             plot3([nf(1),0],[nf(2),0],[nf(3),0],'--','color','b');hold on;
%             plot3([ni(1),0],[ni(2),0],[ni(3),0],'--','color','g');hold on;
            
            %Force vector;
            plot3([fv(1)+xm,xm],[fv(2)+ym,ym],[fv(3)+zm,zm],'--','color','b');hold on;
            
            %manoeuvre
            scatter3(xm,ym,zm,'black');hold on;
            
            axis vis3d
            axis equal
            xlabel('x[m]')
            ylabel('y[m]')
            zlabel('z[m]')
            
            %a=plane_change(30,67,0,10000,700+1736,0,67).run().cplot();
        end
    end
end

%name=plane_change(10,20,5,10000,1900,20,30).run();