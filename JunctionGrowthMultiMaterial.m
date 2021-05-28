
% Material Properties - Constituent #1
den = 2700;         %  density:                     kg/m^3      
hc = 896;           %  heat capacity:               J/kg*K      
vhc = den*hc;       %  volumetric heat capcity:     J/m^3*K
tcb = 167;          %  thermal conductivity:        W/m*K       
tdb = tcb/vhc;      %  thermal diffusivity
mp = 925;           %  melting point:                 K            
emod =70e9;         %  elastic modulus                Pa
vp = 0.33;          %  Poisson ratio                  (-)
em = emod/(1-(vp^2)); % reduced elasctic modulus      Pa

% Material Properties - Consittuent #2
den2 = 8890;         %  density:                     kg/m^3      
hc2 = 385;           %  heat capacity:               J/kg*K      
vhc2 = den2*hc2;       %  volumetric heat capcity:     J/m^3*K
tcb2 = 388;          %  thermal conductivity:        W/m*K       
tdb2 = tcb2/vhc2;      %  thermal diffusivity
mp2 = 1356;           %  melting point:                 K            
emod2 =122.5e9;         %  elastic modulus                Pa
vp2 = 0.33;          %  Poisson ratio                  (-)
em2 = emod2/(1-(vp^2)); % reduced elasctic modulus      Pa

% Material Properties - Sonotrode
tcsono = 6.8;               %  thermal conductivity:        W/m*k      
vhcsono = 2.34e6;           %  volumeteric heat capcity     J/m^3*K     
tdsono = tcsono/vhcsono;    %  thermal diffusivity:    m^2/s
dsono = 0.096;              %  sonotrode diameter           m
sonowidth = 0.0254;         % sonotrode width               m

% Process Parameters
tempamb = 338;          % ambient temperature:          K
amp = 18e-6;            % oscillation amplitude:        m
freq = 20e3;            % oscillation frequency:         s^-1       
strainr = 1e4;          % strain rate:                 s^-1
hpf = tcb/(tcb+sqrt(tdb/tdsono)*tcsono); % heat partition factor
power = 2600;           % avg power consumed by machine                     W
sonovel = 0.0741;       % sonotrode velocity used in experimental work  m/s
normalload = 4500;         % normal load used in experimental work       N 
CoF = power/(sqrt(2)*pi*normalload*freq*amp);      % coefficient of friction is calculated based on the power and process parameters 

% Johnson-Cook Parameters - Consituent #1
A = 236;        
B = 430;        
D = 0.024;      
m = 1.34;
srr = 1;

% Johnson Cook Parameters - Consistuent #2
ACu = 99;
BCu = 253;
DCu = 0.023;
mCu = 0.9;

% Data Spacing
loadstep = 100;                  % number of normal load intervals   increase this to smooth out contours at cost of computation time
loadlimit = 9000;               % maximum normal load (N)
velocitystep = 100;              % number of sonotrode velocity intervals      increase this to smooth out contours at cost of computation time
velocitylimit = 0.1;            % maximum sonotrode velocity (m/s)
xsteps = 10;                     % number of steps in position vector along interface

% Matrix Initialization
temp =zeros(1,xsteps);      
tempmat = zeros(loadstep,velocitystep);
hardmat = zeros(loadstep,velocitystep);
hardmat2 = zeros(loadstep,velocitystep);
presmat = zeros(loadstep,velocitystep);
bmat = zeros(loadstep,velocitystep);
lv = zeros(loadstep,1);                 % vector with range of normal loads
vv = zeros(1,velocitystep);             % vector with range of sonotrode velocities

for load = loadlimit/loadstep : loadlimit/loadstep : loadlimit          % loops through range of loads
    lv(load*loadstep/loadlimit) = load/1000;                            % populates normal load range vector
   
    for v = velocitylimit/velocitystep : velocitylimit/velocitystep : velocitylimit     % loops through sonotrode velocities
    
        b = sqrt( (2*load*dsono)/(pi*sonowidth) *  (((1-vp^2)/emod) + ((1-vp2^2)/emod2))   );         % half width of sonotrode contact region
        xspace = linspace(-b,b,xsteps);                     % position along interface (x) vector
        k = (v/(2*tdb));                                        % defined constant
        pre = (2*sqrt(2)*freq*CoF*amp*hpf*load)/(b*pi*tcb*sonowidth);       % prefactor to temperature integral (contains heat flux)
        f = @(X,x) pre*real(sqrt(1 - (X/b).^2)).*exp(-k.*(x-X)).*besselk(0,k.*abs(x-X));  % function to calculate temperature f(x)
        g = @(x) integral(@(X) f(X,x), -b,b);                   % integrate f over x to get temperature 
            for xi = 1:xsteps                               % loop through x to get temperature at each location
                temp(xi) = tempamb + g(xspace(xi));
            end
           
            presmat(int16(load*loadstep/loadlimit),int16(v*velocitystep/velocitylimit)) =  (2*load)/(pi*b*sonowidth);    % max pressure matrix
            hardmat(int16(load*loadstep/loadlimit),int16(v*velocitystep/velocitylimit))= 3*1e6*(A+B)*(1+(D*log(strainr/srr))) ...  
                *(1-((max(temp)-tempamb)/(mp-tempamb))^m);        % hardness from Johnson-Cook matrix
            hardmat2(int16(load*loadstep/loadlimit),int16(v*velocitystep/velocitylimit))= 3*1e6*(ACu+BCu)*(1+(DCu*log(strainr/srr))) ...  
                *(1-((max(temp)-tempamb)/(mp-tempamb))^mCu);        % hardness from Johnson-Cook matrix
            tempmat(int16(load*loadstep/loadlimit),int16(v*velocitystep/velocitylimit)) = max(temp);    % max temperature matrix     
            vv(int16(v*velocitystep/velocitylimit)) = v*1000;    % populates velocity range vector
    end
    prog = num2str((100*load)/loadlimit);     
    progress = " " + prog + " " + '%' + " "
    
end

    % Contant Area Matrices Calculation
    amat = presmat./hardmat; 
    amat(amat<0) = 1;
    amat(amat > 1) = 1;
    amat2 = presmat./hardmat2; 
    amat2(amat2<0) = 1;
    amat2(amat2 > 1) = 1;
    
    
    areacon = [0.1, 0.15, 0.3, 0.89, 0.999];        % contact area contours shown in output plot
    tempcon = [450, 500, 700,mp, mp2];              % temperature contours
    outputmatrix = zeros(1000,8+2*length(areacon)); % matrix the contains load vs velocity data for contact area and temperature contours
    
[c,w] = contour(vv,lv,amat,areacon);      % creates contour map of contact area matrix
k=1;     % contour line number
col=1;   % index of column containing contour level and number of points
while col<size(c,2) % while less than total columns in c
   s(k).level = c(1,col); 
   s(k).numel = c(2,col); 
   idx=col+1:col+c(2,col);
   s(k).xdata = c(1,idx).'; 
   s(k).ydata = c(2,idx).'; 
   s(k).isopen = abs(diff(c(1,idx([1 end]))))>1e-12 || ...
                 abs(diff(c(2,idx([1 end]))))>1e-12; 
   k=k+1;
   col=col+c(2,col)+1;
end

for k = 1:length(areacon)
    
    xs = linspace(s(k).xdata(1), s(k).xdata(end), 1000);      % populates matrix "output matrix" with area contours x-y data
    y = interp1(s(k).xdata, s(k).ydata, xs);
    outputmatrix(:,(2*k)-1) = xs;
    outputmatrix(:,2*k) = y;
end


[c,w] = contour(vv,lv,amat2,areacon);      % creates contour map of contact area matrix
k=1;     % contour line number
col=1;   % index of column containing contour level and number of points
while col<size(c,2) % while less than total columns in c
   s(k).level = c(1,col); 
   s(k).numel = c(2,col); 
   idx=col+1:col+c(2,col);
   s(k).xdata = c(1,idx).'; 
   s(k).ydata = c(2,idx).'; 
   s(k).isopen = abs(diff(c(1,idx([1 end]))))>1e-12 || ...
                 abs(diff(c(2,idx([1 end]))))>1e-12; 
   k=k+1;
   col=col+c(2,col)+1;
end
    xs = linspace(s(5).xdata(1), s(5).xdata(end), 1000);    
    y = interp1(s(5).xdata, s(5).ydata, xs);
    outputmatrix(:,13) = xs;
    outputmatrix(:,14) = y;
    
    xs = linspace(s(4).xdata(1), s(4).xdata(end), 1000);    
    y = interp1(s(4).xdata, s(4).ydata, xs);
    outputmatrix(:,17) = xs;
    outputmatrix(:,18) = y;

[c,w] = contour(vv,lv,tempmat,tempcon);  
k=1;     % contour line number
col=1;   % index of column containing contour level and number of points
while col<size(c,2) 
   s(k).level = c(1,col); 
   s(k).numel = c(2,col); 
   idx=col+1:col+c(2,col);
   s(k).xdata = c(1,idx).'; 
   s(k).ydata = c(2,idx).'; 
   s(k).isopen = abs(diff(c(1,idx([1 end]))))>1e-12 || ...
                 abs(diff(c(2,idx([1 end]))))>1e-12; 
   k=k+1;
   col=col+c(2,col)+1;
end    
    xs = linspace(s(4).xdata(1), s(4).xdata(end), 1000);    % extracts melting point contour and converts to x-y data in 
    y = interp1(s(4).xdata, s(4).ydata, xs);
    outputmatrix(:,11) = xs;
    outputmatrix(:,12) = y;   
    
    xs = linspace(s(5).xdata(1), s(5).xdata(end), 1000);    % extracts melting point contour and converts to x-y data in 
    y = interp1(s(5).xdata, s(5).ydata, xs);
    outputmatrix(:,15) = xs;
    outputmatrix(:,16) = y;  
    
close all

% Plot calculations - red contours show melting point of each material,
% blue contours show paerameters that lead to complete bonding
figure
plt = plot(...
outputmatrix(:,9),outputmatrix(:,10),'b',outputmatrix(:,11),outputmatrix(:,12),'red'...
,outputmatrix(:,13),outputmatrix(:,14),'b',outputmatrix(:,15),outputmatrix(:,16),'red');
set(plt,'LineWidth',2);   
set(gca,'linewidth',1.5);
set(gca, 'FontSize', 16);
xlabel('sonotrode velocity (mm/s)', 'FontSize',20);
ylabel('normal load (kN)', 'FontSize',20);
text(max(outputmatrix(:,11)),max(outputmatrix(:,12)),'melt boundary #1','Color','red','FontSize',14);
text(max(outputmatrix(:,15))-4,max(outputmatrix(:,16))-0.5,'melt boundary #2','Color','red','FontSize',14);
text(max(outputmatrix(:,9)),max(outputmatrix(:,10)),'Ar/An-#1 =1','FontSize',14);
text(max(outputmatrix(:,13)),max(outputmatrix(:,14)),'Ar/An-#2 = 1','FontSize',14);




