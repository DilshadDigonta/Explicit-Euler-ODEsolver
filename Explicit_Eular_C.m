function Explicit_Eular_C()

close all
clear all

s = pi/18;
% Solving the ODE for
% 0<t<4 with x1(0)=s and x2(0)=0 [range of t and x ]

%Own Euler Solver
[t1,x1]=MyEulerSys(@dxdt,[0 4],[s 0],0.2);                                   %[range of t,x and h]
[t2,x2]=MyEulerSys(@dxdt,[0 4],[s 0],0.1);
[t3,x3]=MyEulerSys(@dxdt,[0 4],[s 0],0.05);

%Matlab Solver
[tmat,xmat]=ode23(@dxdt,[0 4],[s 0]);
plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'g-',tmat,xmat(:,1),'ro-')
legend('h=0.2','h=0.1','h=0.05','ode23()');                                 %[ploting the result]
xlabel('t');
ylabel('x');
title("Exp Euler");


function [t,x]=MyEulerSys(ODEfunc,tspan,x0,h)
%This function uses Euler method to solve a
%system of ODEs
t=tspan(1):h:tspan(2);
x(1,:)=x0;
for n=1:length(t)-1
    x(n+1,:)=x(n,:)+ODEfunc(t(n),x(n,:))'*h;
end
 
 
function xp=dxdt(t,x)                                                       %[function for the ODE, canonical form] 
xp(1)=x(2);
xp(2)=-16.35.*x(1);
xp=xp';