% Finite Difference Method - ELEC 4700 Assignment 2

close all
clear 
format short
clc

global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

% Constants
m0 = 9.11e-31; %kg
me = 0.26*m0;
q = 1.602e-19; %C
kB = 1.38066e-23; %J/K
nm = 1e-9; %nanometre
ps = 1e-12; %picosecond

n=100; % nodes
l=1.5; 
w=1.0; 

x=linspace(0,l,n); %length
y=linspace(0,w,n); % width
V=ones(n)*0.1;
Vo = 0.1;

steps = 100;
% limit = 1e-4;
% error = 1;

% Boundary conditions 1a:
V(n,1:n)=0; % Bottom
V(1,1:n)=0; % Top
V(1:n,1)=Vo; % left
V(1:n,n)=0; % right

for k = 1:steps
% while error > limit
%     
%     Vprev = V;
       
    for i=2:n-1 % skipping BC in iterations
        for j=2:n-1
            V(i,j)=(V(i+1,j)+V(i-1,j)+V(i,j-1)+V(i,j+1))/4;
        end
    end  
    
%     error=max(max(abs(Vprev-V)));
    
    figure(1)
    subplot(2,1,1),pcolor(x,y,V),shading interp, colormap;
    title('2D Voltage Plots')
    xlabel('x')
    ylabel('y')
       
    subplot(2,1,2),
    [Ex,Ey] = gradient(V);
    quiver(-Ey,-Ex,10);
    title('Electric Field (V/m)')
    
%     pause(0.01)
end

% Boundary conditions 1b:
V(n,1:n)=0; % Bottom
V(1,1:n)=0; % Top
V(1:n,1)=Vo; % left
V(1:n,n)=Vo; % right

for k = 1:steps
       
    for i=2:n-1 % skipping BC in iterations
        for j=2:n-1
            V(i,j)=(V(i+1,j)+V(i-1,j)+V(i,j-1)+V(i,j+1))/4;
        end
    end       
    
    figure(2)
    subplot(2,1,1),pcolor(x,y,V),shading interp, colormap;
    title('2D Voltage Plots')
    xlabel('x')
    ylabel('y')
    colorbar;
       
    subplot(2,1,2),
    [Ex,Ey] = gradient(V);
    quiver(-Ey,-Ex,10);
    title('Current Density')
    
    pause(0.01)
end

figure(3)
surf(x,y,V,'EdgeColor','none')     
xlabel('x'),ylabel('y'),zlabel('V(x,y)')
title('Numerical Solution Mesh')
rotate3d on

a = l*n;
b = w*n/2;
add = 0;
       
% for i=1:n 
%     for j=1:n           
%         for k = 1:10
% 
%             add = add +(1/k)*cosh(k*pi*i/a)/cosh(k*pi*b/a)*sin(k*pi*j/a);
%             V(i,j)=4*Vo/pi*add;
% 
%             figure(4)
%             surf(x,y,V,'EdgeColor','none')     
%             xlabel('x'),ylabel('y'),zlabel('V(x,y)')
%             title('Analytical Solution Mesh')
% 
%             pause(0.01)
%         end 
%     end
% end       

% Boundary conditions 2:
V(n,1:n)=0; % Bottom
V(1,1:n)=0; % Top
V(1:n,1)=Vo; % left
V(1:n,n)=0; % right

%2nd set of BCs
Boxes{1}.X = [0.5 1.0]; %top box
Boxes{1}.Y = [0.7 1.0];
%         Boxes{1}.BC = 0.0;

Boxes{2}.X = [0.5 1.0]; %bottom box
Boxes{2}.Y = [0.0 0.3];
%         Boxes{2}.BC = 0.0;

%defining conduction matrix
sigma_in = 1e-2;
sigma_out = 1;
sig = ones(n);
ibc = i*l/n;
jbc = j*w/n;

% defining conduction in boxes
for i = floor(Boxes{1}.X(1,1)*n/l):ceil(Boxes{1}.X(1,2)*n/l)
    for j = 1:Boxes{2}.Y(1,2)*n
        sig(i,j) = sigma_in;
    end
    for j = Boxes{1}.Y(1,1)*n:n
        sig(i,j) = sigma_in;
    end
end


for k = 1:steps
% while error > limit
%     
%     Vprev = V;
       
    for i=2:n-1 % skipping BC in iterations
        for j=2:n-1
            V(i,j)=sig(i,j)*(V(i+1,j)+V(i-1,j)+V(i,j-1)+V(i,j+1))/4;            
        end
    end  
    
%     error=max(max(abs(Vprev-V)));
    
    figure(5)
    subplot(2,1,1),pcolor(x,y,V),shading interp, colormap;
    title('2D Voltage Plots')
    xlabel('x')
    ylabel('y')
    hold on

    %box lines
    bx1 = [0.5 0.5];
    bx2 = [1 1];
    bx3 = [0.5 1];
    by1 = [0.7 1];
    by2 = [0 0.3];
    by3 = [0.3 0.3];
    by4 = [0.7 0.7];

    %box plot
    plot(bx1,by1,'k')
    plot(bx1,by2,'k')
    plot(bx2,by1,'k')
    plot(bx2,by2,'k')
    plot(bx3,by3,'k')
    plot(bx3,by4,'k')
    hold off

    subplot(2,1,2),
    [Ex,Ey] = gradient(V);
    quiver(-Ey',-Ex',50);
    title('Electric Field (V/m)')
    hold on 
    
    %box lines
    bx1 = [0.5 0.5];
    bx2 = [1 1];
    bx3 = [0.5 1];
    by1 = [0.7 1];
    by2 = [0 0.3];
    by3 = [0.3 0.3];
    by4 = [0.7 0.7];

    %box plot
    plot(bx1,by1,'k')
    plot(bx1,by2,'k')
    plot(bx2,by1,'k')
    plot(bx2,by2,'k')
    plot(bx3,by3,'k')
    plot(bx3,by4,'k')
    hold off
    
    if mod(k,steps/10)==0
        pause(0.05)
    end
end

figure(6) 
subplot(2,1,1),pcolor(sig'),
title('Conductivity Plot')
xlabel('x'),ylabel('y')

subplot(2,1,2),pcolor(V)
title('Electric Field')
xlabel('x'),ylabel('y')

figure(7)
Area = 500:100:2000;
I = 100./Area;
plot(Area,I)
title('Bottle Neck Area vs Current')
xlabel('Bottle Neck Area'),ylabel('Current')

figure(8)
sig = 00:0.1:1;
I = 100*sig;
plot(sig,I)
title('Conductivity vs Current')
xlabel('Conudctivity'),ylabel('Current')

