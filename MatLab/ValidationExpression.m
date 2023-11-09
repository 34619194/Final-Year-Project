sizeR = 'small';
shape = 'butterfly';

N=1000;
phi = linspace(0,2*pi,N);
fung = zeros(1,N);
vphi = zeros(1,N);

for i=1:N
    [fung(i),~,~] = getGs(phi(i), sizeR, shape);
    vphi(i) = mod(fung(i),2*pi);
end
varphi = umwrap(vphi);

function [delta,Ddelta, DDdelta,DDDdelta] = getDELTAs(phi,shape)end
function [alpha,Dalpha,DDalpha] = getALPHAs(phi,shape)end
function [g,Dg,DDg] = getGs(phi,sizeR,shape)end

