function [delta,Ddelta,DDdelta,DDDdelta] = getDELTAs(b,c,phi)
%getDELTAs
%    [DELTA,Ddelta,DDdelta,DDDdelta] = getDELTAs(B,C,PHI)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    05-Nov-2023 12:31:12

t2 = phi.*2.0;
t3 = cos(t2);
delta = b-c.*t3;
if nargout > 1
    t4 = sin(t2);
    Ddelta = c.*t4.*2.0;
end
if nargout > 2
    DDdelta = c.*t3.*4.0;
end
if nargout > 3
    DDDdelta = c.*t4.*-8.0;
end
