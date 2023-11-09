
%define object class
CR = CircleRobot;
%We Choose which program to run:
program = 0;
if (program == 0)
%Simulate alpha, beta , gamma

to_  = 0; %Initail time 
tEnd_ = 10; %final time
x0_ = [0; 1];
[t_,x_] = CR.simABG(t0_ , tEnd_ , x0_);

elseif (program ==1)
%Simulate the full circle system

t0 = 0;
tEnd = 5;
x0 = [0; pi/2; 0; 0] ;
[t ,x] = CR.sim(t0, tEnd, x0);

end 