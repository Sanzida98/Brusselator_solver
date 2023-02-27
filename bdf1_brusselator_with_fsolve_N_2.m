
%memory cleanup:
clear all;
close all;
clc;

%parameters:
N=2; %Spatial points
%initial condition:
j=1;

%for u in the odd rows:
for i=1:2:2*N 
    if j<=N
         y0(i,1) = 1+sin((2*pi/(N+1))*j);
    end
    j=j+1;
end

%for v in the even rows:
for i=2:2:N+2
    y0(i,1)=3;
end

