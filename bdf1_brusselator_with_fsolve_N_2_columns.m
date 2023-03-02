
%memory cleanup:
clc;
clear all;
close all;
clearvars;

%parameters:
N=2; %Spatial points
a=1/50; %function constant

%initial condition:
j=1;
for i=1:2:2*N 
    %Odd rows:
    if j<=N
       y0(i,1) = 1+sin((2*pi/(N+1))*j);% initial condition
    end
    j=j+1;
    %Even rows:
    y0(i+1,1)=3;
end

%time:
t0=0;
tfinal=2;
h=0.1; %step size
t=t0:h:tfinal; %time vector
n=length(t); %time vector points

%Solution vector:
y=zeros(2*N,n); %expands in column with time
y(:,1)=y0; %initiating solution vector

%fsolve:
options = optimset('Display','off','Jacobian','on'); %fsolve optimization
for i=2:length(t)
    y0(:,i) = y(:,i-1); %previous step solution
    y(:,i) = fsolve( @(y) BEJ(y,a,N, y0(:,i), t(i),h,@ode,@ode_j), y0(:,i), options);
end

%Brusselator system function for N=2:
function f=ode(y,a,N)
   f1=1+ (y(1,:).^2 .*y(2,:))-(4 .*y(1,:))+((a*(N+1)^2).*(1-2.*y(1,:)+y(3,:)));
   g1=3.*y(1,:)-(y(1,:).^2 .*y(2,:))+((a*(N+1)^2).*(3-2.*y(2,:)+y(4,:)));
   f2=1+(y(3,:).^2 .*y(4,:))-(4.*y(3,:))+((a*(N+1)^2).*(y(1,:)-2 .*y(3,:)+1));
   g2=(3.*y(3,:))-(y(3,:).^2 .*y(4,:))+((a*(N + 1)^2).*(y(2,:)-2 .*y(4,:)+3));
   f=[f1;g1;f2;g2]; %system ode
end

%System Jacobian for N=2:
function j=ode_j(y,a,N)
     j1=[(2.*y(1,:).*y(2,:))-4-(2*a*(N+1)^2), (y(1,:).^2), a*(N+1)^2, 0];
     j2=[(3-2.*y(1,:).*y(2,:)), (-y(1,:).^2 -2*a*(N+1)^2), 0, a*(N+1)^2];
     j3=[(a*(N+1)^2), 0, (2.*y(3,:).*y(4,:)-4-2*a*(N+1)^2), (y(3,:).^2)];
     j4=[0, (a*(N+1)^2), (3-2.*y(3,:).*y(4,:)), (-y(3,:).^2-2*a*(N+1)^2)];
     j=[j1;j2;j3;j4];
end

%bdf1 solver function:
function [F,J] = BEJ(y,a,N,y0,t,h,ode,ode_j)
       %F=zeros(2*N,1);
       F=y-y0-h*ode(y,a,N); %function input for fsolve in terms of bdf1
       if nargout>1 %check if Jacobian is supplied to fsolve
       J=eye(2*N)-h*ode_j(y,a,N); %function jacobian
       end
end