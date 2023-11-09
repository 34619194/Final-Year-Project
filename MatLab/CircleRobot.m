classdef CircleRobot
%DynamicsofCircleRobot
properties
m=3e-3;%massoftheball
Jf=1.581e-3;%Momentifinteruiaoftheframe
Jb=5.48e-7;%Momentofinteriaoftheball
g=[0;9.81;0];%Gravity
R_b=16.55e-3;%Radiusoftheball
r_f=12.5e-3;%Distancebetweenplates
Rf=0.1;%Radiusofframeinmeters
R;%Effectiveradius
Rd;%Radiusoforigiotoballscentre
c=0.49;%ScalarforVHC
B_an=[0 1];%Annihilatormatrix
phiApprox;%Approxiamtionofphi
end

methods
%Constructor
function obj = CircleRobot()
obj.R= sqrt(obj.R_b^2 - obj.r_f^2);
obj.Rd=obj.R+obj.Rf;
N=100;
phi=linspace(0,2*pi,N);
g=zeros(1,N);
for i=1:N
[g(i),~,~]=getGs(obj,phi(i));
end
varphi=mod(g,2*pi);
obj.phiApprox=spline(varphi,phi);
end

%Alpha,Beta,Gamma
function[t,x]=simABG(obj,t0,tEnd,x0)
options=odeset('RelTol',1e-5,'AbsTol',1e-6);
tspan=linspace(t0,tEnd,1000);
[t,x]=ode23(@(t,y)ABG_EOM(obj,y),tspan,x0,options);
end

function dx=ABG_EOM(obj,x)
vphi=x(1);
Dvphi=x(2);
[alpha,beta,gamma]=getABG(obj,vphi);
DDvphi=-(beta*Dvphi^2 + gamma)/alpha;
dx=[Dvphi;DDvphi];
end

function[alpha,beta,gamma]=getABG(obj,varphi)
[Phi,DPhi,DDPhi]=getPhis(obj,varphi);
Dq=[Phi;DPhi];
[M,C,G]=getMCG(obj,Dq);
alpha=obj.B_an*M*DPhi;
beta=obj.B_an*(C*DPhi+M*DDPhi);
gamma=obj.B_an*G;
end

function[Phi,DPhi,DDPhi]=getPhis(obj,varphi)
[Theta,DTheta,DDTheta]=getVHC(obj,varphi);
Phi=[Theta;varphi];
DPhi=[DTheta; 1];
DDPhi=[DDTheta; 0];
end

function[Theta,DTheta,DDTheta]=getVHC(obj,varphi)
Theta=varphi-obj.c*sin(2*varphi);
DTheta=1-2*obj.c*cos(2*varphi);
DDTheta=4*obj.c*sin(2*varphi);
end

function[t,x]=simCR(obj,t0,tEnd,x0)
options=odeset('RelTol',1e-5,'AbsTol',1e-6);
tSteps=linspace(t0,tEnd,100);
[t,x]=ode23(@(t,y)Butterfly_EOM(obj,y),tSteps,x0,options);
end

function dx=Butterfly_EOM(obj,x)
Dq=x(3:4);
[M,C,G]=getMCG(obj,x);
DDq=M\(-C*Dq-G);
dx=[Dq;DDq];
end

%sharedfunctions
function[M,C,G]=getMCG(obj,x)
%parameters
q1=x(1);q2=x(2);Dql=x(3);Dq2=x(4);
Jf=obj.Jf;g=obj.g;
Jb=obj.Jb;m=obj.m;
R=obj.R;Rd=obj.Rd;
[z,Dz,Pi,DPi]=getVariables(obj,q1,q2);

%MakeM
m11=Jf+Jb+m*Rd^(2);
m12=-(m*Rd^(2)+(Jb/R)*Rd);
m21=m12;
m22=(m+(Jb/R^(2)))*Rd^(2);
M=[m11,m12;m21,m22];

%MakeC
c11=0;
c12=0;
c21=0;
c22=0;
C=[c11,c12;c21,c22];

%MakeG
g1=dot(m*g,DPi*Rd*z);
g2=dot(m*g,Pi*Rd*Dz);
G=[g1;g2];
end

function[z,Dz,Pi,DPi]=getVariables(obj,theta,varphi)
phi=ppval(obj.phiApprox,varphi);

Pi=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
DPi=[-sin(theta) -cos(theta) 0; cos(theta) -sin(theta) 0; 0 0 0];

z=[sin(phi);cos(phi);0];
Dz=[cos(phi);-sin(phi);0];
end

 %Symbolic Functions 
    function [delta ,Ddelta ,DDdelta, DDDdelta] = getDELTAs(obj , phi)
        a = obj.a; b = obj.b;
        t2 = phi.*2.0;
        t3 = cos(t2);
        delta = a-b.*t3;
        t4 = sin(t2);
        Ddelta = b.*t4.*2.0;
        DDdelta = b.*t3.*4.0;
        DDDdelta = b.*t4.*-8.0;
        if strcmpi(obj.shape , 'circle')
            delta = obj.Rf;
            Ddelta = 0;
            DDdelta = 0;
            DDDdelta = 0;
        end 
    end
function [alpha ,Dalpha ,DDalpha] = getALPHAs(obj , phi)
        [delta ,Ddelta ,DDdelta ,DDDdelta] = getDELTAs(obj ,phi);
        t2 = delta;
        t3 = sin(phi);
        t4 = Ddelta;
        t5 = cos(phi);
        alpha = atan((t2.*t3 - t4.*t5)./(t2.*t5 + t3.*t4));
        if(phi > pi/2 && phi <3*pi/2)
            alpha = alpha + pi;
        end
        if(phi >3*pi/2)
            alpha = alpha +2*pi;
        end 
        t6 = t4.^2;
        t7 = t2.^2;
        t8 = t6 + t7;
        t9 = DDdelta;
        Dalpha = (t6.*2.0 + t7 - t2.*t9)./t8;
        t10 = DDDdelta;
        DDalpha = -1.0./t8.^2.*(t2.*t4.*t9.^2.*-2.0 + t2.*t4.*t6.*2.0 + t2.*t6.*t10 + t2.*t7.*t10 + t4.*t6.*t9 - t4.*t7.*t9.*3.0);
    end
    function [g, Dg, DDg] = getGs(obj , phi)
        R = obj.R;
        [delta ,Ddelta, DDdelta ,~] = getDELTAs(obj ,phi);
        [alpha ,Dalpha , DDalpha] = getALPHAs(obj , phi);
        t2 = delta;
        t3 = alpha;
        g = atan2((R.*sin(t3)+t2.*sin(phi)) ,(R.*cos(t3) + t2.*cos(phi)));

        if(nargout >1)
            t4 = t2.^2;
            t5 = R.^2;
            t6 = phi-t3;
            t7 = cos(t6);
            t8 = Dalpha;
            t9 = R.*t2.*t7.*2.0;
            t10 = t4+t5+t9;
            t11 = 1.0./t10;
            t12 = Ddelta;
            t13 = sin(t6);
            t14 = DDalpha;
            t15 = t2.*t12.*2.0;
            t16 = R.*t7.*t12.*2.0;
            t17 = t5.*t8;
            t18 = R.*t12.*t13;
            t19 = R.*t2.*t7.*t8;
            t20 = R.*t2.*t7.*t8;
            t21 = t4+t17+t18+t19+t20;
            Dg = t11.*t21;
            DDg = t11.*(t15 + t16 + t5.*t14 - R.*t2.*t13 + R.*t13.*DDdelta  + R.*t2.*t7.*t14 + R.*t2.*t8.^2.*t13)-1.0./t10.^2.*t21.*(t15+t16+R.*t2.*t13.*(t8-1.0).*2.0);
        end 
    end
end
end


