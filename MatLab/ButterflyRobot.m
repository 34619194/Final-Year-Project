classdef ButterflyRobot
%Dynamics of Butterflyrobot
properties
    m       = 3e-3;         %mass of the ball
    Jf      = 1.581e-3;     %Moment if interuia of the frame
    Jb      = 5.48e-7;      %Moment of interia of the ball
    g       = [0; 9.81; 0];   %Gravity 
    R_b     = 16.55e-3;     %Radius of the ball
    r_f     = 12.5e-3;      %Distance between plates
    Rf      = 0.1;          %Radius of frame in meters 
    R;                      %Effective radius 
    k_hat   = [0;0;1];      %z directional vector
    B       = [1;0];        %Coupling matrix 
    B_an    = [0 1];        %Annihilator matrix 
    a       = 0.1095;       % constant for butterfly frame 
    b       = 0.0405;       % scalar for butterfly frame 
    c       = 0.49;         % Scalar for VHC
    shape = 'butterfly';    % choose shape form 
    phiApprox;              %Approxiamtion of phi
    I = eye(3);             %3x3 Identity matrix 
    Q = [0 1 0;             %Matrix for z' = Qz
        -1 0 0;
        0 0 0]; 
end

methods 
    % Constructor 
    function obj = ButterflyRobot()
        obj.R = sqrt(obj.R_b^2 - obj.r_f^2);
        N = 100;
        phi = linspace(0.2*pi,N);
        g = zeros(1,N);
        for i = 1:N
            g(i) = getGs(obj , phi(i));
        end
        vphi = mod(g, 2*pi);
        varphi = unwrap(vphi);
        obj.phiApprox = spline(varphi, phi);
    end

    % Alpha, Beta, Gamma 
    function [t,x] = simABG(obj , t0 , tEnd , x0)
        options = odeset('RelTol' , 1e-5, 'AbsTol',1e-6);
        tspan = linspace(t0 , tEnd , 500);
        [t,x] = ode23(@(t,y)ABG_EOM(obj,y) , tspan, x0, options);
    end
    function dx = ABG_EOM(obj , x)
        vphi =x(1);
        Dvphi = x(2);
        [alpha, beta, gamma] = getABG(obj,vphi);
        DDvphi = -(beta*Dvphi^2 + gamma)/alpha;
        dx = [Dvphi;DDvphi];
    end

    %Butterfly Robot 
    function [t,x,Fn,Fs,u] = simBR(obj,t0,tEnd,x0)
        options = odeset('RelTol' , 1e-5, 'AbsTol',1e-6);
        tSteps = linspace(t0,tEnd,500);
        [t,x] = ode23(@(t,y)Butterfly_EOM(obj,y) , tSteps, x0, options);
        n = length(t);
        u = zeros(n,1);
        for i=1:n
            u(i) = getActuatorInput(obj ,x(i,:)');
        end
        Fn = zeros(n,1);
        Fs = zeros(n,1);
        for i=1:n
            [Fn(i) ,Fs(i)] = getForces(obj ,x(i,:)');
        end
    end

    function dx = Butterfly_EOM(obj ,x)
        [M,C,G] = getMCG(obj,x);
        u = getActuatorInput(obj,x);
        Dq = x(3:4);
        DDq = M\(-C*Dq - G + obj.B*u);
        dx = [Dq;DDq];
    end

    function u = getActuatorInput(obj,x)
        %Genralized coordinates 
        theta = x(1); Dtheta = x(3);
        varphi = x(2); Dvarphi = x(4);

        %VHC
        Theta = varphi - obj.c*sin(2*varphi);
        DTheta = 1 - 2*obj.c*cos(2*varphi);
        DDTheta = 4*obj.c*sin(2*varphi);

        %Dynamics 
        [M,C,G] = getMCG(obj,x);

        %reduced dynamics 
        [alpha , beta , gamma] = getABG(obj , varphi);
        
        %Acceleration
        DDvarphi = -(beta/alpha)*Dvarphi^2 -(gamma/alpha);

        %Accutation
        u = M(1,:)*[DDTheta*Dvarphi^2 + DTheta*DDvarphi; DDvarphi] + C(1,:)*[DTheta*Dvarphi;Dvarphi] + G(1);

    end

    function [Fn,Fs] = getForces(obj,x)
        Dx = Butterfly_EOM(obj, x);
        Dtheta = Dx(1); DDtheta =Dx(3);
        Dvarphi = Dx(2); DDvarphi = Dx(4);

        x_mod = [x(1) ; x(2) ; 0; 0];
        [~,rhoVec , tau ,n , kappa ,Ds ,DDs ,Pi, ~] = getVariables(obj,x_mod);
     
        %Normal force
        Fn = obj.m*(Pi*n)'*(2*Pi*n*Ds*Dtheta*Dvarphi + Pi*kappa*Ds^(2)*Dvarphi^(2)+cross(DDtheta*obj.k_hat , Pi*rhoVec) + obj.g);
        
        %Friction force 

        Fs = obj.m*(Pi*tau)'*(Pi*kappa*Ds^(2)*Dvarphi^(2)+Pi*tau*DDs*Dvarphi^(2) + Pi*tau*Ds*DDvarphi+cross(DDtheta*obj.k_hat , Pi*rhoVec) + obj.g);

    end 

    %shared functions 
    function[M,C,G] = getMCG(obj , x)
        %parameters 
        q1 = x(1) ; q2 = x(2) ; Dq1 = x(3) ; Dq2 = x(4) ;
        m = obj.m ; Jf = obj.Jf ; Jb = obj.Jb ; g = obj.g ; R = obj.R ; 
        k_hat = obj.k_hat ; B = obj.B ; a = obj.a ; b = obj.b;
        [rho ,rhoVec , tau , ~ , kappa , Ds ,DDs, Pi ,DPi] = getVariables (obj , x );
        
        %Make M
        m11 = Jf + Jb + m*rho^(2);
        m12 = (dot(m*k_hat , cross(rhoVec , tau))-(Jb/R))*Ds;
        m21 = m12;
        m22 = (m + (Jb/R^(2)))*Ds^(2);
        M = [m11,m12;m21,m22];
       
        %Make C
        c11 = dot(m*Ds*rhoVec , tau)*Dq2;
        c12 = dot(m*Ds*rhoVec , tau)*Dq1 + ((dot(m*k_hat , cross(rhoVec , tau)) - (Jb/R))*DDs + dot(m*k_hat , cross(rhoVec , kappa))*Ds^(2))*Dq2;
        c21 = -dot(m*Ds*rhoVec , tau)*Dq1;
        c22 = (m + (Jb/R^(2)))*Ds*DDs*Dq2;
        C = [c11,c12;c21,c22];


        %Make G
        g1 = dot(m*g , DPi*rhoVec);
        g2 = dot(m*g , Pi*tau*Ds);
        G = [g1;g2];
    end 
    
    function [rho ,rhoVec ,tau ,n ,kappa ,Ds ,DDs ,Pi ,DPi] = getVariables(obj , x)
        
        %Genral coordinates and angles
        theta  = x(1);
        varphi = x(2);
        phi = ppval(obj.phiApprox, varphi);
        [alpha, Dalpha ,DDalpha] = getALPHAs(obj , phi);

        %Vectors needed for calculations 
        [deltaVec ,DdeltaVec ,DDdeltaVec] = getDeltaVecs(obj ,phi);

        %Tangens vectors
       tau = DdeltaVec/norm(DdeltaVec);

       %Normal Vectors
       n = [sin(alpha) ;cos(alpha) ;0];
       Dn = [cos(alpha)*Dalpha ;- sin(alpha)*Dalpha ;0];
       DDn = [-sin(alpha)*Dalpha^2 + cos(alpha)*DDalpha ;-cos(alpha)*Dalpha^2 - sin(alpha)*DDalpha ; 0];

       %Curvature Vector 

       kappa_f = norm(cross(DdeltaVec ,DDdeltaVec))/((norm(DdeltaVec))^3);
       kappa_const = kappa_f/(1-obj.R*kappa_f);
       kappa = kappa_const*n;

       %vector from origin to intersection point: ball/frame
       rhoVec = deltaVec + obj.R*n;
       rho = norm(rhoVec);

       %Arc length, differianted once and twice
       [~, Dg, DDg] = getGs(obj , phi);
       h = DdeltaVec + obj.R*Dn;
       Dh = DDdeltaVec + obj.R*DDn;
       Ds = (norm(h))*(1/Dg);
       DDs = (((dot(h,Dh)/norm(h))*(1/Dg)) - ((norm(h))*(DDg/Dg^2)))*(1/Dg);

       Pi = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
       DPi = [-sin(theta) -cos(theta) 0; cos(theta) -sin(theta) 0; 0 0 0];
    end 

    function [deltaVec ,DdeltaVec ,DDdeltaVec] = getDeltaVecs(obj , phi)
        [delta ,Ddelta ,DDdelta ] = getDELTAs(obj , phi);
        z = [sin(phi) ;cos(phi) ; 0];
        deltaVec = delta*z;
        DdeltaVec = (Ddelta*obj.I + delta*obj.Q)*z;
        DDdeltaVec = (DDdelta*obj.I +2*Ddelta*obj.Q - delta*obj.I)*z;
    end

    function [alpha ,beta ,gamma] = getABG(obj , varphi)
        [Phi ,DPhi ,DDPhi] = getPhis(obj ,varphi);
        Dq = [Phi ; DPhi];
        [M ,C ,G ] = getMCG(obj ,Dq);
        alpha = obj.B_an*M*DPhi;
        beta = obj.B_an*(C*DPhi + M*DDPhi);
        gamma = obj.B_an*G;
    end 

    function [Phi ,DPhi ,DDPhi] = getPhis(obj , varphi)
        [Theta ,DTheta ,DDTheta ] = getVHC(obj , varphi);
        Phi = [Theta; varphi ];
        DPhi = [DTheta; 1 ];
        DDPhi = [DDTheta; 0];
    end 

    function [Theta ,DTheta ,DDTheta] = getVHC( obj , varphi)
        Theta = varphi - obj.c*sin(2*varphi);
        DTheta = 1 - 2*obj.c*cos(2*varphi);
        DDTheta = 4*obj.c*sin(2*varphi);      
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






