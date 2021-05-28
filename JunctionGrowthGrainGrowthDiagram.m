clear all

% Material Properties - Nanocrystalline Constituent
ae1 = 1.264;         %  activation energy           eV          
gge1 = 3;            %  grain growth exponent       –           
mob1 = 1.30259E-11;  %  grain boundary mobility     m^gge1/s        
d01 = 15e-9;         %  initial grain size          nm
    
den1 = 8900;         %  density                     kg/m^3      
hc1 = 446;           %  heat capacity               J/kg*K      
vhc1 = den1*hc1;     %  volumetric heat capcity     J/m^3*K
tc1 = 82.9;          %  thermal conductivity        W/m*K       
td1 = tc1/vhc1;      %  thermal diffusivity         m/s^2
mp1 = 1725;          %  melting point               K          
emod1 = 200e9;       %  elastic modulus             Pa
vp1 = 0.31;          %  Poisson's ratio             –
    

% Material Properties -  Constituent #2
den2 = 2710;         %  density                     kg/m^3     
hc2 = 904;           %  heat capacity               J/kg*K     
vhc2 = den2*hc2;     %  volumetric heat capcity     J/m^3*K
tc2 = 222;           %  thermal conductivity        W/m*K
td2 = tc2/vhc2;      %  thermal diffusivity         m/s^2
mp2 = 930;           %  melting point               K                   
emod2 = 69e9;        %  elastic modulus             Pa
vp2 = 0.33;          %  Poisson's ratio             –


% Properties - Sonotrode 
tcs = 25.5;          %  thermal conductivity        W/m*k    
dens = 8080;         %  density                     kg/m^3  
hcs = 481;           %  heat capcity                J/kg*K     
vhcs = dens*hcs;     %  volumetric heat capcity     J/m^3*K
tds = tcs/vhcs;      %  thermal diffusivity         m^2/s
ds = 0.0973582;      %  sonotrode diameter          m
sw = 0.0254;         %  sonotrode width             m
                      

% Material Properties – Base Plate   
tcbp = 60.2;         %  thermal conductivity        W/m*k  
hcbp = 481;          %  heat capcity                J/kg*K  
denbp = 7872;        %  density                     kg/m^3 
vhcbp = denbp*hcbp;  %  volumetric heat capcity     J/m^3*K
tdbp = tcbp/vhcbp;   %  thermal diffusivity         m^2/s

% Material Properties - Driver foils
tcD = 9.35;          %  thermal conductivity        W/m*k     
vhcD = 2.34e6;       %  volumeteric heat capcity    J/m^3*K     

% Johnson-Cook Parameters for Softer Constituent 
A1 = 83;
B1 = 76;
C1 = 0.46;
D1 = 0.49;
m1 = 0.51;
p1 = 0.1;
srr1 = 1;

% Process Parameters / Misc. 
tamb = 330;          %  ambient temperature         K
amp = (32.5e-6)/2;   %  oscillation amplitude       m
freq = 20e3;         %  oscillation frequency       Hz      
sr = 1e4;            %  strain rate                 s^-1
cof = 0.2;           %  coefficient of friction     –  
kb = 8.62e-5;        %  Boltzmann's constant        eV/K


% Grid Data Spacing - higher step values increase resolution of contours
% but increase run time
lstep = 100;        % # of intervals in normal load space (–)
llim = 15000;       % Maximum normal load (N)
vstep = 100;        % # of intervals in normal load space (–)
vlim = 0.1;         % Maximum sonotrode velocity (m/s)

% Matrix/Vector Initialization
xsteps = 200;           % # of position intervals in travel direction of sonotrode (–)
temp =zeros(1,xsteps);          
tempmat = zeros(lstep,vstep);
grainmat = zeros(lstep,vstep);
hardmat = zeros(lstep,vstep);
presmat = zeros(lstep,vstep);
lv = zeros(lstep,1);
vv = zeros(1,vstep);

% Cycles through all normal loads and velocities
for load = llim/lstep : llim/lstep : llim
    lv(load*lstep/llim) = load;             % Creates vector of all load values in grid
    for v = vlim/vstep : vlim/vstep : vlim
    vv(int16(v*vstep/vlim)) = v;            % Creates vector of all velocity values in grid
        b =  sqrt( (2*load*ds)/(pi*sw)*(((1-vp2^2)/emod2) + ((1-vp1^2)/emod1)) );  % Calculate Herztian contact width
        
        d1 = 134e-6;    % Nanocrystalline Constituent Thickness (m)
        d2 = 110e-6;    % Constituent #2 Thickness (m)
        dD = 150e-6;    % Driver Foil Thickness (m)
        
        timedif = (2*b/v);          % Sonotrode contact time (s)
        dbp = sqrt(tdbp*timedif);   % Heat diffusion distance into build plate (m)

        % Weighting Factors
        f2 = d2/(d2+d1+dD+dbp); 
        f1 = d1/(d2+d1+dD+dbp);
        fD = dD/(d2+d1+dD+dbp);
        fbp = dbp/(d2+d1+dD+dbp);
        
        % Effective Thermal Properties
        tceffb = ((f2/tc2)+(f1/tc1)+(fD/tcD)+(fbp/tcbp))^-1;
        vhceffb = (f2*vhc2)+(f1*vhc1)+(fD*vhcD)+(fbp*vhcbp);
        tdeffb = tceffb/vhceffb;
       
        hpf = 0.5; % Heat partition factor (–)
      
        xspace = linspace(-5*b,2*b,xsteps);                     % Position vector for distance along travel direction       
        k = (v/(2*tdeffb));                                     % Convienent grouping of parameters
        pre = (2*sqrt(2)*freq*cof*amp*hpf*load)/(pi*tceffb*sw); % Convienent grouping of parameters
        
        % Defining Temperature Functions
        f = @(X,x) (pre/b)*real(sqrt(1 - (X/b).^2)).*exp(-k.*(x-X)).*besselk(0,k.*sqrt((x-X).^2));
        g = @(x) integral(@(X) f(X,x), -b,b); 
            for xi = 1:xsteps
                temp(xi) = tamb + g(xspace(xi));
            end
           
            time = (xspace./v)-xspace(1)/v;                 % Convert position vector to time vector
            di = mob1*(1./temp).*exp(-ae1./(kb.*temp));     % Integrand of grain growth equation
            z = cumtrapz(time,di);                          % Numerically integrate grain growth equation
            grain = ( (d01)^gge1 + z).^(1/gge1);            % Calculate grain size as function of time vector
            
            % Input calculated data to pressure, hardness, temperature, and
            % grain size matrices
            presmat(int16(load*lstep/llim),int16(v*vstep/vlim)) =  (2*load)/(pi*b*sw);
            hardmat(int16(load*lstep/llim),int16(v*vstep/vlim)) = (3*1e6*(A1+B1)*(C1+(D1*log(sr/srr1))^p1) ...  
                *(1-((max(temp)-tamb)/(mp2-tamb))^m1));
            tempmat(int16(load*lstep/llim),int16(v*vstep/vlim)) = max(temp);
            grainmat(int16(load*lstep/llim),int16(v*vstep/vlim)) = max(grain)/1e-9;
    end
    progress = load*lstep/llim
end
close all

% Contant Area Matrix Calculation
amat = presmat./hardmat; 
amat(amat<0) = 1;
amat(amat > 1) = 1;

% Define Contour Values
areacon = [0.8, 1];
tempcon = [500, mp2];
graincon = [20,25,35,50,100];

OutputMatrix = zeros(1000,2*length(graincon)+4);

% Extract Ar/An=1 Contour from Contact Area Matrix
[d,w] = contour(vv,lv,amat,areacon); 
k=1;     % contour line number
col=1;   % index of column containing contour level and number of points
while col<size(d,2) % while less than total columns in c
   s(k).level = d(1,col); 
   s(k).numel = d(2,col); 
   idx=col+1:col+d(2,col);
   s(k).xdata = d(1,idx).'; 
   s(k).ydata = d(2,idx).'; 
   s(k).isopen = abs(diff(d(1,idx([1 end]))))>1e-12 || ...
                 abs(diff(d(2,idx([1 end]))))>1e-12; 
   k=k+1;
   col=col+d(2,col)+1;
end
    xw1 = linspace(s(2).xdata(1), s(2).xdata(end),1000);
    [s(2).xdata, index] = unique(s(2).xdata);
    
    % Store v vs. L data for Ar/An = 1 Contour in Output Matrix
    OutputMatrix(:,2) = interp1(s(2).xdata, s(2).ydata(index), xw1);
    OutputMatrix(:,1) = xw1;

    
% Extract Grain Contour from Contact Area Matrix
[c,w] = contour(vv,lv,grainmat,graincon);  
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
   

for k = 1:length(graincon)
    x = linspace(s(k).xdata(1), s(k).xdata(end), 1000);
    y = interp1(s(k).xdata, s(k).ydata, x);
    OutputMatrix(:,3+ 2*k) = x;
    OutputMatrix(:,4+ 2*k) = y;
end


   
% Extract Melting Point Contou from Contact Area Matrix
[a,w] = contour(vv,lv,tempmat,tempcon); 
k=1;     % contour line number
col=1;   % index of column containing contour level and number of points
while col<size(a,2) % while less than total columns in c
   s(k).level = a(1,col); 
   s(k).numel = a(2,col); 
   idx=col+1:col+a(2,col);
   s(k).xdata = a(1,idx).'; 
   s(k).ydata = a(2,idx).'; 
   s(k).isopen = abs(diff(a(1,idx([1 end]))))>1e-12 || ...
                 abs(diff(a(2,idx([1 end]))))>1e-12; 
   k=k+1;
   col=col+a(2,col)+1;
end
    OutputMatrix(:,3) = linspace(s(2).xdata(1), s(2).xdata(end), 1000);
    OutputMatrix(:,4) = interp1(s(2).xdata, s(2).ydata, OutputMatrix(:,3));
    
for m = 1:length(graincon)  
    for m1 = 1:1000
    
        if OutputMatrix(m1,(2*m)+4) >= OutputMatrix(m1,4)
              OutputMatrix(m1,(2*m)+4) = NaN;  
        end
    end
end
    
close all

% Plot Ar/An=1 contour, metling point, and grain size contours 
figure
plt = plot(OutputMatrix(:,1),OutputMatrix(:,2),'blue', OutputMatrix(:,3),OutputMatrix(:,4),'red',OutputMatrix(:,5),...
    OutputMatrix(:,6),'green',OutputMatrix(:,7),OutputMatrix(:,8),'green',OutputMatrix(:,9),...
    OutputMatrix(:,10),'green',OutputMatrix(:,11),OutputMatrix(:,12),'green',OutputMatrix(:,13),OutputMatrix(:,14),'green');
set(plt,'LineWidth',2);    
set(gca,'linewidth',1.5)
set(gca, 'FontSize', 16)
xlabel('sonotrode velocity', 'FontSize',20)
ylabel('weld pressure (MPa)', 'FontSize',20)




