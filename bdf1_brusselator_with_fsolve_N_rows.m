%bdf1+fsolve solver for Brusselator system of N:

%memory cleanup:
clc;
clear all;
close all;
clearvars;

%System parameters:
N=input('Type the value of N:   '); %Spatial points
a=1/50; %function constant

%initial condition:
j=1;
for i=1:2:2*N 
    %Odd rows:
    if j<=N
       y0(1,i) = 1+sin((2*pi/(N+1))*j);% initial condition
    end
    j=j+1;
    %Even rows:
    y0(1,i+1)=3;
end

%time array:
t0=0;
tfinal=2;
h=0.1; %step size
t1=t0:h:tfinal; %time vector
n=length(t1); %time vector points

%Solution vector:
y=zeros(n,2*N); %expands in column with time
y(1,:)=y0; %initiating solution vector

%fsolve:
options = optimset('Display','off','Jacobian','on'); %fsolve optimization
for i=2:length(t1)
    y0(i,:) = y(i-1,:); %previous step solution
    y(i,:) = fsolve( @(y) BEJ(y,a,N, y0(i,:), t1(i),h,@ode,@ode_j), y0(i,:), options);
end

%Load data from brusscode.m:
load("plot_solution.mat"); %solution vector
load("plot_space.mat"); %space array
load("plot_time.mat"); %time array

%plot of solver solution:
u1 = y(:,1:2:end);
x1 = (1:N)/(N+1); %defining x
figure;
surf(x1,t1,u1,'FaceColor','g');
hold on;

%plot of brusscode.m solution:
surf(x,t,u,'FaceColor','b','Marker','o');
colormap('cool');

function f=ode(y,a,N)
%for first two rows: 
i=1;
f(1,i)=1+ (y(:,i).^2 *y(:,i+1))-4*y(:,i)+(a*(N+1)^2).*(1-2*y(1,i)+y(i+2));
f(1,i+1)=3*y(1,i)- (y(1,i).^2 .*y(i+1))+(a*(N+1)^2)*(3-2*y(1,i+1)+y(i+3));

%for last two rows:
i=2*N-1;
f(1,i)=1+(y(:,i).^2 .*y(:,i+1))-4.*y(:,i)+(a*(N+1)^2)*(y(:,i-2)-2*y(:,i)+1);
f(1,i+1)=3.*y(:,i)-y(:,i).^2 .*y(:,i+1)+(a*(N+1)^2)*(y(:,i-1)-2*y(:,i+1)+3);

%for rest of the rows in between:
j=1;
if N>2
     for i=2:N-1
         if j<2*N
             j=j+2;
             f(1,j)=1+y(:,i)^2 *y(:,i+1)-4*y(:,i)+(a*(N+1)^2)*(y(:,i-2)-2*y(:,i)+y(:,i+2));
             f(1,j+1)=3*y(i)-y(i)^2 *y(i)+c*(y(i-1)-2*y(i)+y(i+1)); 
         end
     end
end
end

%System Jacobian for N=2:
function j=ode_j(y,a,N)
     j1=[(2.*y(:,1).*y(:,2))-4-(2*a*(N+1)^2), (y(:,1).^2), a*(N+1)^2, 0];
     j2=[(3-2.*y(:,1).*y(:,2)), (-y(:,1).^2 -2*a*(N+1)^2), 0, a*(N+1)^2];
     j3=[(a*(N+1)^2), 0, (2.*y(:,3).*y(:,4)-4-2*a*(N+1)^2), (y(:,3).^2)];
     j4=[0, (a*(N+1)^2), (3-2.*y(:,3).*y(:,4)), (-y(:,3).^2-2*a*(N+1)^2)];
     j=[j1;j2;j3;j4];
end

%bdf1 solver function:
function [F,J] = BEJ(y,a,N,y0,t1,h,ode,ode_j)
       %F=zeros(2*N,1);
       F=y-y0-h*ode(y,a,N); %function input for fsolve in terms of bdf1
       if nargout>1 %check if Jacobian is supplied to fsolve
       J=eye(2*N)-h*ode_j(y,a,N); %function jacobian
       end
end

