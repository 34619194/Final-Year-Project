m = 3e-3;
Jf = 1.581e-3;
Jb = 5.48e-7;
g =[0;9.81;0];
R_b = 16.55e-3;
r_f = 12.5e-3;
R = sqrt(R_b^2 - r_f^2);
k_hat = [0;0;1];
B = [1;0];
B_an = [0 1];
c = 0.49;

N=1000;
phi = linspace(0,2*pi,N);
fung = zeros(1,N);
varphi = zeros(1,N);
VHC = zeros(1,N);
VHC_check = zeros(1,N);
VHC_DTheta = zeros(1,N);
VHC_DDTheta = zeros(1,N);
Alpha = zeros(1,N);
Beta = zeros(1,N);
Gamma = zeros(1,N);
DGamma = zeros(1,N);

for i=1:N

    [fung(i),Dg,DDg] = getGs(phi(i));
    varphi(i) = mod(fung(i),2*pi);

    [alpha,Dalpha,DDalpha] = getALPHAs(phi(i));
    [deltaVec,DdeltaVec,DDdeltaVec] = getDeltaVecs(phi(i));
    tau = DdeltaVec/norm(DdeltaVec);
    n = [sin(alpha); cos(alpha); 0];
    Dn = [cos(alpha)*Dalpha; -sin(alpha)*Dalpha; 0];
    DDn = [-sin(alpha)*Dalpha^2 + cos(alpha)*DDalpha; -cos(alpha)*Dalpha^2 - sin(alpha)*DDalpha; 0];
    rho = deltaVec + R*n;
    h = DdeltaVec + R*n;
    Dh = DDdeltaVec + R*DDn;
    Ds = (norm(h))*(1/Dg);
    DDs = (((dot(h,Dh)/norm(h))*(1/Dg))-((norm(h))*(DDg/Dg^2)))*(1/Dg);

    Theta = varphi(i) - c*sin(2*varphi(i));
    DTheta = 1 - 2*c*cos(2*varphi(i));
    DDTheta = 4*c*sin(2*varphi(i));
    VHC(i) = Theta;
    VHC_check(i) = (-Ds*(m + Jb/R^2))/(dot(m*k_hat,cros(rho,tau)) - Jb/R);

    Pi = [cos(Theta) - sin(Theta) 0; sin(Theta) cos(Theta) 0; 0 0 1];
    DPi = [-sin(Theta) -cos(Theta) 0; cos(Theta) -sin(Theta) 0; 0 0 0];

    Alpha(i) = Ds*(dot(m*k_hat,cross(rho,tau)) - Jb/R)*DTheta + Ds^(2)*(m + Jb/R^2);
    Beta(i) = Ds*(dot(m*k_hat,cross(rho,tau))-Jb/R)*DDTheta - dot(m*Ds'*tau,rho*DTheta.^(2)) + (m + Jb/R^2)*Ds*Dss;
    Gamma(i) = dot(m*Ds*g,Pi*tau);
    DGamma(i) = m*DDs*dot(g,Pi*tau) + m*Ds*dot(g,(DPi*DTheta*tau + Pi*((DDdeltaVec/norm(DdeltaVec))-((DdeltaVec*dot(DdeltaVec,DDdeltaVec))/(norm(DDdeltaVec)^3)))));
end

asymCheckOne = Gamma./Alpha;
asymCheckTwo = Beta./Alpha;
[deltaVecEq1,DdeltaVecEq1,DDdeltaVecEq1] = getDeltaVecs(0);
[deltaVecEq2,DdeltaVecEq2,DDdeltaVecEq2] = getDeltaVecs(pi/2);
[alphaEq1,DalphaEq1,~] = getALPHAs(0);
[alphaEq2,DalphaEq2,~] = getALPHAs(pi/2);
[~,DgEq1,~]= getGs(0);
[~,DgEq2,~]= getGs(pi/2);
rhoEqOne = deltaVecEq1 + R*[sin(alphaEq1); cos(alphaEq1); 0];
rhoEqTwo = deltaVecEq2 + R*[sin(alphaEq2); cos(alphaEq2); 0];
tauEqOne = DdeltaVecEq1/norm(DdeltaVecEq1);
tauEqTwo = DdeltaVecEq2/norm(DdeltaVecEq2);
DsEqnOne = (norm(DdeltaVecEq1 + R*[cos(alphaEq1)*DalphaEq1; -sin(alphaEq1)*DalphaEq1; 0]))*(1/DgEq1);
DsEqnTwo = (norm(DdeltaVecEq2 + R*[cos(alphaEq2)*DalphaEq2; -sin(alphaEq2)*DalphaEq2; 0]))*(1/DgEq2);

aOptions = (-2:2);
for i = aOptions
    DThetaEqOne = 1 - 2*i*cos(2*0);
    DThetaEqTwo = 1 - 2*i*cos(2*(pi/2));
    AlphaEqPointOne(i+3) = DsEqnOne.*(dot(m*k_hat,cross(rhoEqOne,tauEqOne)) - Jb/R).*DThetaEqOne + DsEqnOne.^(2).*(m + Jb/R^2);
    AlphaEqPointTwo(i+3) = DsEqnTwo.*(dot(m*k_hat,cross(rhoEqTwo,tauEqTwo)) - Jb/R).*DThetaEqTwo + DsEqnTwo.^(2).*(m + Jb/R^2);
end 

AlphaNull = zeros(length(aOptions));
P1 = InterX([aOptions;AlphaEqPointOne],[aOptions,AlphaEqPointTwo]);
P2 = InterX([aOptions;AlphaNull],[aOptions,AlphaEqPointOne]);
P3 = InterX([aOptions;AlphaNull],[aOptions,AlphaEqPointTwo]);
Px = [p1(1);P2(1);P3(1)];
Py = [p1(2);P2(2);P3(2)];

lambdaSignDecider = DGamma./Alpha;

function [deltaVec,DdeltaVec,DDdeltaVec] = getDeltaVecs(phi)
[delta,Ddelta,DDdelta,~] = getDELTAs(phi);
z = [sin(phi); cos(phi); 0];
Q = [0 1 0; -1 0 0; 0 0 0];
I = eye(3);
deltaVec = delta*z;
DdeltaVec = (Ddelta*I + delta*Q)*z;
DDdeltaVec = (DDdelta*I + 2*Ddelta*Q - delta*I)*z;
end 

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





