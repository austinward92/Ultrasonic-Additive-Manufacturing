
% Material Properties - Feedstock
den = 2700;         %  density:                     kg/m^3      
hc = 896;           %  heat capacity:               J/kg*K      
vhc = den*hc;       %  volumetric heat capcity:     J/m^3*K
tcb = 167;          %  thermal conductivity:        W/m*K      
tdb = tcb/vhc;      %  thermal diffusivity
mp = 925;           %  melting point:               K          

% Material Properties - Sonotrode
tcs = 6.8;          %  thermal conductivity:        W/m*k    
vhcs = 2.34e6;      %  volumeteric heat capcity     J/m^3*K     
tds = tcs/vhcs;     %  thermal diffusivity:         m^2/s
ds = 0.096;         %  sonotrode diameter           m
sw = 0.0254;        %  sonotrode width              m
emod =70e9;         %  elastic modulus              Pa
vp = 0.33;          %  Poisson's ratio              –
em = emod/(1-(vp^2));  %  Reduced elastic modulus   Pa

% Process Parameters
tamb = 293;         % ambient temperature:          K
amp = 28e-6;        % oscillation amplitude:        m
freq = 20e3;        % oscillation frequency:        s^-1
cof = 0.2;          % coefficient of friction:      -
sr = 1e4;           % strain rate:                  s^-1

hpf = tcb/(tcb+sqrt(tdb/tds)*tcs);  % Heat partition factor

% Johnson-Cook Parameters - Feedstock
A = 236;
B = 430;
D = 0.024;
m = 1.34;
srr = 1;

% Grid Data Spacing - higher step values increase resolution of contours
% but increase run time
lstep = 100;     % # of intervals in normal load space (–)
llim = 10000;   % Maximum normal load (N)
vstep = 100;     % # of intervals in normal load space (–)
vlim = 0.1;     % Maximum sonotrode velocity (m/s)

xsteps = 10;     % # of position intervals in travel direction of sonotrode (–)

% Matrix/Vector Initialization
temp =zeros(1,xsteps);
tempmat = zeros(lstep,vstep);
hardmat = zeros(lstep,vstep);
presmat = zeros(lstep,vstep);
bmat = zeros(lstep,vstep);
lv = zeros(lstep,1);
vv = zeros(1,vstep);

% Cycles through all normal loads and velocities
for load = llim/lstep : llim/lstep : llim
    lv(load*lstep/llim) = load/1000;          % Creates vector of all load values in grid
    for v = vlim/vstep : vlim/vstep : vlim
    vv(int16(v*vstep/vlim)) = v*1000;         % Creates vector of all velocity values in grid
        
        b = sqrt((2*load*ds)/(pi*sw*em));       % Calculate Herztian contact width
        xspace = linspace(-b,b,xsteps);          % Position vector for distance along travel direction  
        k = (v/(2*tdb));                            % Convienent grouping of parameters
        pre = (2*sqrt(2)*freq*cof*amp*hpf*load)/(b*pi*tcb*sw);   % Convienent grouping of parameters
        
        % Defining Temperature Functions
        f = @(X,x) pre*real(sqrt(1 - (X/b).^2)).*exp(-k.*(x-X)).*besselk(0,k.*abs(x-X));
        g = @(x) integral(@(X) f(X,x), -b,b);
            for xi = 1:xsteps
                temp(xi) = tamb + g(xspace(xi));
            end
           
            % Populate output matrices
            presmat(int16(load*lstep/llim),int16(v*vstep/vlim)) =  (2*load)/(pi*b*sw);
            hardmat(int16(load*lstep/llim),int16(v*vstep/vlim))= 3*1e6*(A+B)*(1+(D*log(sr/srr)))*(1-((max(temp)-tamb)/(mp-tamb))^m);
            tempmat(int16(load*lstep/llim),int16(v*vstep/vlim)) = max(temp);
    end
    progress = load*lstep/llim
end

% Contant Area Matrix Calculation
    amat = presmat./hardmat; 
    amat(amat<0) = 1;
    amat(amat > 1) = 1;

% Define Contour Values
    areacon = [0.1, 0.2, 0.3,0.5,0.999];
    tempcon = [450, 500, 700,mp];
    
[c,w] = contour(vv,lv,amat,areacon);  
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
   
OutputMatrix = zeros(1000,2+2*length(areacon));

for k = 1:length(areacon)
    
    xs = linspace(s(k).xdata(1), s(k).xdata(end), 1000);
    y = interp1(s(k).xdata, s(k).ydata, xs);
    OutputMatrix(:,(2*k)-1) = xs;
    OutputMatrix(:,2*k) = y;
end

[c,w] = contour(vv,lv,tempmat,tempcon);  
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

    xs = linspace(s(4).xdata(1), s(4).xdata(end), 1000);
    y = interp1(s(4).xdata, s(4).ydata, xs);
    OutputMatrix(:,11) = xs;
    OutputMatrix(:,12) = y;


close all

% green is Ar/An=1 (complete bonding), red is melting boundary, blue
% contours are rest of Ar/An values
figure
plt = plot(OutputMatrix(:,1),OutputMatrix(:,2),'b', OutputMatrix(:,3),OutputMatrix(:,4),'b',OutputMatrix(:,5),OutputMatrix(:,6),'b',OutputMatrix(:,7),OutputMatrix(:,8),'b',OutputMatrix(:,9),OutputMatrix(:,10),'g',OutputMatrix(:,11),OutputMatrix(:,12),'red');
set(plt,'LineWidth',2);    
set(gca,'linewidth',1.5)
set(gca, 'FontSize', 16)
xlabel('sonotrode velocity (mm/s)', 'FontSize',20)
ylabel('normal load (kN)', 'FontSize',20)


