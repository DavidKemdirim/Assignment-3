% MC/FD - ELEC 4700 Assignment 1

% Monte Carlo Modelling of Electron Transport - ELEC 4700 Assignment 1

close all
clear 
format short
clc

% Constants
m0 = 9.11e-31; %kg
me = 0.26*m0;
q = 1.602e-19; %C
kB = 1.38066e-23; %J/K
nm = 1e-9; %nanometre
ps = 1e-12; %picosecond

% Dimensions
% Elec = 1; % simulates for 1 particle
Elec = 5; % simulates for # particles at once
xdim = 200; %nm
ydim = 100; %nm and need to make same length
x = zeros(Elec,1)*nm;
y = zeros(Elec,1)*nm; %need to make same length 
% Reg = zeros(xdim*1e9,ydim*1e9)*nm;  % semiconductor region

% Specifics
InitialTemp = 300; % in [K]
vtFirst = sqrt(2*kB*InitialTemp/me); %1.8701e5 in m/s
% Tmn = 0.2*ps; 
Tmn = 0.2*ps/3;

% Simulation
steps = 10; 
% steps = 1000;
% t = zeros(1,steps);
% Temp = zeros(1,steps);
dt = nm/vtFirst; %5.347e-15s
Temp = zeros(Elec,1);
TempNew = zeros(Elec,1);
vtNew = zeros(Elec,1);
t = zeros(1,Elec);
dz = zeros(Elec,1); % straight path
dx = zeros(Elec,1); % new x dif
dy = zeros(Elec,1); % new y dif 

on = 1;
off = 0;
scatter = on; %scattering on/off switch
% scatter = off; %scattering on/off switch
collide = 0; % counter for collisons
p = 1 - exp(-dt/Tmn); % probability for scatter threshold
re = 0; % default rebound velocity factor
% re = -0.5; % Bounces back at half its initial speed

xbc = 0;
xCompare = 0;
ybc = 0;
yCompare = 0;

for i = 1:steps    
   
    if i == 1
        t(i,:) = dt;
        Temp(:,i) = InitialTemp;
        vt = vtFirst;        
        
        Boxes = {};
        
        %Test BCs
%         Boxes{1}.X = [0.9 1.1]*100*nm; %top box
%         Boxes{1}.Y = [0.8 1.0]*100*nm;
% %         Boxes{1}.BC = 0.0;
%         
%         Boxes{2}.X = [0.9 1.1]*100*nm; %bottom box
%         Boxes{2}.Y = [0.0 0.2]*100*nm;
% %         Boxes{2}.BC = 0.0; 

        %2nd set of BCs
        Boxes{1}.X = [0.8 1.2]*100*nm; %top box
        Boxes{1}.Y = [0.6 1.0]*100*nm;
%         Boxes{1}.BC = 0.0;
        
        Boxes{2}.X = [0.8 1.2]*100*nm; %bottom box
        Boxes{2}.Y = [0.0 0.4]*100*nm;
%         Boxes{2}.BC = 0.0;
        
        e=1;
        while e <= Elec %random positions for each particle
                        
            randx = round(rand*xdim)*nm;
            randy = round(rand*ydim)*nm;
            
            % check to see if in box
            if Boxes{1}.X(1,1)<randx && randx<Boxes{1}.X(1,2) &&...
               Boxes{1}.Y(1,1)<randy && randy<Boxes{1}.Y(1,2) ||...%firstbox
               Boxes{2}.X(1,1)<randx && randx<Boxes{2}.X(1,2) &&...
               Boxes{2}.Y(1,1)<randy && randy<Boxes{2}.Y(1,2) %second box
                                
            else
                x(e,i) = randx;
                y(e,i) = randy;
                e=e+1;
            end
                      
        end 
        
    else
        t(i,:) = t(i-1,:) + dt;
        r = rand(1,Elec); % Chance of scattering
        
        % Particle direction
        if i == 2 % initial angle, regardless of scattering
            
            for e = 1:Elec
                theta = rand*360; % random angle in deg
                dz(e,i) = vtFirst*Tmn; % straight path
                dx(e,i) = sind(theta)*dz(e,i); % new x dif
                dy(e,i) = cosd(theta)*dz(e,i); % new y dif 
                x(e,i) = x(e,i-1) + dx(e,i);
                y(e,i) = y(e,i-1) + dy(e,i); 
                Temp(e,i) = Temp(e,i-1);
                vtNew(e,i) = vtFirst;
            end

        else
            
            if scatter == on % yes scatter
            
                for e = 1:Elec 

                    if r(1,e) < p   % check each electron for scatter %

    %                     scale = sqrt(kB*Temp(e,i-1)/me); %scaling factor
    %                     dof = 3; %degrees of freedom
                        vtNew(e,i) = vt+randn*vt/sqrt(2); % vt in place of scale
    %                     vtNew(e,i) = scale*chi2rnd(dof);
    %                     vtNew(e,i) = vtFirst*(rand/100+0.995);

                        theta = rand*360; % random angle in deg
                        dz(e,i) = vtNew(e,i)*Tmn; % straight path
                        dx(e,i) = sind(theta)*dz(e,i); % new x difference
                        dy(e,i) = cosd(theta)*dz(e,i); % new y difference 

                        % updating position, velocity, temperature
                        TempNew(e,i) = me*vtNew(e,i).^2/(2*kB); 
                        Temp(e,i) = TempNew(e,i); 
                        x(e,i) = x(e,i-1) + dx(e,i); 
                        y(e,i) = y(e,i-1) + dy(e,i); 
    %                     vtNew = vt;

                        collide = collide + 1; 

                    else % no scattering, just continue along inital path 
                        dx(e,i) = x(e,2)-x(e,1);
                        dy(e,i) = y(e,2)-y(e,1); 
                        dz(e,i) = sqrt(dx(e,i)^2 + dy(e,i)^2);
                        x(e,i) = x(e,i-1) + dx(e,i); 
                        y(e,i) = y(e,i-1) + dy(e,i);
    %                     vtNew = vt;
                        Temp(e,i) = Temp(e,i-1);
                    end

                end

            else % no scattering, just continue along inital path

                for e = 1:Elec
                    dx(e,i) = x(e,2)-x(e,1);
                    dy(e,i) = y(e,2)-y(e,1); 
                    dz(e,i) = sqrt(dx(e,i)^2 + dy(e,i)^2);
                    x(e,i) = x(e,i-1) + dx(e,i); 
                    y(e,i) = y(e,i-1) + dy(e,i);
                end
        %             vtNew = vt;
                    Temp(e,i) = Temp(e,i-1);

            end
                
            % Boundary Conditions
            for e = 1:Elec

                if x(e,i) < 0 || x(e,i) > xdim*nm % loop on left/right

                    xbc = xbc + 1;
    %                     theta = rand*360; % random angle in deg
    %                     dz(e,i) = vtNew(e,i)*Tmn; % straight path
    %                     dx(e,i) = sind(theta)*dz(e,i); % new x difference
    %                     dy(e,i) = cosd(theta)*dz(e,i); % new y difference 

                    if x(e,i) < 0
                        x(e,i) = x(e,i) + 2e-7;    
    %                     xNew(e,i) = x(e,i) + 2e-7;
                    elseif x(e,i) > xdim*nm
                        x(e,i) = x(e,i) - 2e-7;
    %                     xNew(e,i) = x(e,i) + 2e-7;
                    end

                end

                if y(e,i) < 0 || y(e,i) > ydim*nm % reflect off top/bottom         

                    if y(e,i) < 0 && dy(e,i) < 0 
                        
                        while y(e,i) < 0 
                            y(e,i) = y(e,i) - dy(e,i)/10;
                        end
                        
                    elseif y(e,i) > ydim*nm && dy(e,i) > 0                            
                        
                        while y(e,i) > ydim*nm
                            y(e,i) = y(e,i) - dy(e,i)/10;
                        end
                        
                    end
            
                    dy(e,i) = cosd(-theta)*dz(e,i); % refelcting y dif 
                    y(e,i) = y(e,i-1) + dy(e,i); 

                end


            end
        
        end
                     
                        
    end

    %box lines
    bx1 = [80 80]*nm;
    bx2 = [120 120]*nm;
    bx3 = [80 120]*nm;
    by1 = [60 100]*nm;
    by2 = [0 40]*nm;
    by3 = [40 40]*nm;
    by4 = [60 60]*nm;
    
%     1c-i)  Particle Plot   
    for e = 1:Elec
        figure(1)
        if xbc > xCompare && ybc > yCompare         
            plot(xNew(e,:),yNew(e,:));
            xbc = xCompare;
            ybc = yCompare;
        elseif xbc < xCompare && ybc > yCompare
            plot(x(e,:),yNew(e,:));
            ybc = yCompare;
        elseif xbc > xCompare && ybc < yCompare
            plot(xNew(e,:),y(e,:));
            xbc = xCompare;
        end
             
        plot(x(e,:),y(e,:));
        
        grid on;
        hold on; 
    end
    %box plot
%     plot(bx1,by1,'k')
%     plot(bx1,by2,'k')
%     plot(bx2,by1,'k')
%     plot(bx2,by2,'k')
%     plot(bx3,by3,'k')
%     plot(bx3,by4,'k')

    xlim([0 xdim*nm])
    ylim([0 ydim*nm])
    xlabel('X (m)')
    ylabel('Y (m)')
    title(['Time Passed t: ', num2str(t(i)/ps), ...
        'ps Collsions: ', num2str(collide)]) 

%     % 1c-ii) Temp plot  
%     for e = 1:Elec
%         figure(2)
%         plot(t(:,e),Temp(e,:));
%         grid on;
%         hold on;
%     end 
%     xlabel('Time (s)')
%     ylabel('Temperature (K)')
%     title(['Average Temperature: ', num2str(mean(Temp(:,i))), ...
%         'K Max T: ', num2str(max(Temp,[],'all')),'K Min T: ',...
%         num2str(min(Temp,[],'all')),'K']) % change T        
         
    pause(0.02)
end

display('Seconds Passed', num2str(t(i)));
finalTemp = num2str(mean(Temp(:,i)))

scale = sqrt(kB*InitialTemp/me);
h = scale*chi2rnd(3,1,10000); 
figure(3)
histogram(vtNew,10)
xlabel('Thermal Velocity (m/s)')
ylabel('Count')
title(['Velocity Distribution - Average Velocity: ',...
    num2str(mean([mean(h),median(h)])/2), ' m/s']) 


