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
tfinal=100;
h=0.01; %step size
t1=t0:h:tfinal; %time vector
n=length(t1); %time vector points

%Solution vector:
y=zeros(n,2*N); %expands in column with time
y(1,:)=y0; %initiating solution vector
tic;
%fsolve:
options = optimset('Display','off','Jacobian','on'); %fsolve optimization
for i=2:length(t1)
    y0(i,:) = y(i-1,:); %previous step solution
    y(i,:) = fsolve( @(y) BEJ(y,a,N, y0(i,:), t1(i),h,@ode,@ode_j), y0(i,:), options);
end
toc;
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

%Brusselator system:
function f=ode(y,a,N)
  %for first two rows: 
  f=zeros(1,2*N); %function array
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
     for i=3:2:2*N-2
         if j<2*N-2
             j=j+2;
             f(1,j)=1+(y(:,i).^2 .*y(:,i+1))-4*y(:,i)+(a*(N+1)^2)*(y(:,i-2)-2*y(:,i)+y(:,i+2));
             f(1,j+1)=3*y(:,i)-(y(:,i).^2 .*y(:,i+1))+(a*(N+1)^2)*(y(i-1)-2*y(:,i+1)+y(i+3)); 
         end
     end
  end
end

%System Jacobian for N=2:
function j=ode_j(y,a,N)
     j=zeros(2*N,2*N);

     %first two rows:
     i=1;
     %first row:
     j(i,i)=2.*y(:,i).*y(:,i+1)-4-2*a*(N+1)^2; %df1/du1
     j(i,i+1)=y(:,i).^2; %df1/dv1
     j(i,i+2)=a*(N+1)^2; %df1/du2

     %second row:
     k=2;
     j(k,i)=3-2.*y(:,i).*y(:,i+1); %dg1/du1
     j(k,i+1)=(-y(:,i).^2)-2*a*(N+1)^2; %dg1/dv1
     j(k,i+3)=a*(N+1)^2; %dg1/dv2

     %rest of the rows:
     %odd rows nonzero points:
     %dfi/du(i-1):
     k=1;
     for i=3:2:2*N
         if k<N-2
            j(i,k)=a*(N+1)^2;
            k=k+2;
        end
     end

    %dfi/dui:
    k=1;
    for i=3:2:2*N 
        if k<N
           k=k+2;
           j(i,i)=2*y(:,k).*y(:,k+1)-4-2*a*(N+1)^2;
       end
    end

    %dfi/du(i+1):
    k=3;
    if N>2
       for i=3:2:2*N-2
           if k<N
              k=k+2;
              j(i,k)=a*(N+1)^2;
           end
       end
    end

    %dfi/dvi:
    k=1;
    for i=3:2:2*N
        if k<2*N
           k=k+2;
           j(i,i+1)=y(:,k).^2;
        end
    end

   %Even rows nonzero points:
   %dgi/dui:
   k=1;
   for i=4:2:2*N
      if k<2*N
         k=k+2;
         j(i,k)=3-2*y(:,i-1).*y(:,i);
     end
   end

   %dgi/dv(i-1):
   k=0;
   for i=4:2:2*N
       if k<2*N-2
          k=k+2;
          j(i,k)=a*(N+1)^2;
       end
   end

   %dgi/dvi:
   for i=4:2:2*N
       j(i,i)=(-y(:,i-1).^2)-2*a*(N+1)^2;
   end

   %dgi/dv(i+1):
   k=4;
   if N>4
       for i=4:2:2*N-2
           if k<2*N
              k=k+2;
              j(i,k)=a*(N+1)^2;
           end
       end
   end

end

%bdf1 solver function:
function [F,J] = BEJ(y,a,N,y0,t1,h,ode,ode_j)
       %F=zeros(2*N,1);
       F=y-y0-h*ode(y,a,N); %function input for fsolve in terms of bdf1
       if nargout>1 %check if Jacobian is supplied to fsolve
          J=eye(2*N)-h*ode_j(y,a,N); %function jacobian
       end
end
