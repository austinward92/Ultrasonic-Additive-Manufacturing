clear all

% Material Properties - Constituent #1
den1 = 8900;            %  density                     kg/m^3   
hc1 = 446;              %  heat capacity               J/kg*K   
vhc1 = den1*hc1;        %  volumetric heat capcity     J/m^3*K
tc1 = 82.9;             %  thermal conductivity        W/m*K  
td1 = tc1/vhc1;         %  thermal diffusivity         m/s^2
mp1 = 1725;             %  melting point               K  
emod1 = 200e9;          %  elastic modulus             Pa
vp1 = 0.31;             %  Poisson's ratio             –

% Material Properties -  Constituent #2
den2 = 2710;         %  density                     kg/m^3     
hc2 = 904;           %  heat capacity               J/kg*K     
vhc2 = den2*hc2;     %  volumetric heat capcity     J/m^3*K
tc2 = 222;           %  thermal conductivity        W/m*K
td2 = tc2/vhc2;      %  thermal diffusivity         m/s^2
mp2 = 930;           %  melting point               K                   
emod2 = 69e9;        %  elastic modulus             Pa
vp2 = 0.33;          %  Poisson's ratio             –

% Diffusion constants
DL = 6.5e7;         %   Lattice interdiffusivity prefactor              (um^2/s)   
QL = 2.735;         %   Lattice interdiffusion activation energy        (eV)
DGB = 1e6;          %   Grain boundary interdiffusivity prefactor       (um^2/s)
QGB = 1.16;         %   Grain boundary interdiffusion activation energy (eV)

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

% Process Parameters / Misc. 
tamb = 330;          %  ambient temperature         K
amp = (32.5e-6)/2;   %  oscillation amplitude       m
freq = 20e3;         %  oscillation frequency       Hz      
sr = 1e4;            %  strain rate                 s^-1
cof = 0.2;           %  coefficient of friction     –  
kb = 8.62e-5;        %  Boltzmann's constant        eV/K

% Material Properties - Driver foils
tcD = 9.35;          %  thermal conductivity        W/m*k     
vhcD = 2.34e6;       %  volumeteric heat capcity    J/m^3*K   

% Grid Data Spacing - higher step values increase resolution of contours
% but increase run time
lstep = 40;             % # of intervals in normal load space (–)
llim = 10000;           % Maximum normal load (N)
vstep = 40;             % # of intervals in normal load space (–)
vlim = 0.08;            % Maximum sonotrode velocity (m/s)

% Matrix/Vector Initialization
xsteps = 500;               % # of position intervals in travel direction of sonotrode (–)
temp =zeros(1,xsteps);
dmat = zeros(lstep,vstep);
lv = zeros(lstep,1);
vv = zeros(1,vstep);

% Cycles through all normal loads and velocities
for load = llim/lstep : llim/lstep : llim
    lv(load*lstep/llim) = load;             % Creates vector of all load values in grid
    for v = vlim/vstep : vlim/vstep : vlim
            vv(int16(v*vstep/vlim)) = v;      % Creates vector of all velocity values in grid
            b =  sqrt( (2*load*ds)/(pi*sw)*(((1-vp2^2)/emod2) + ((1-vp1^2)/emod1)) );
        
        d1 = 134e-6;    % Nanocrystalline Constituent Thickness (m)
        d2 = 110e-6;    % Constituent #2 Thickness (m)
        dD = 150e-6;    % Driver Foil Thickness (m)
        
        timedif = (2*b/v);  % Sonotrode contact time (s)
        dbp = sqrt(tdbp*timedif);   % Heat diffusion distance into build plate (m)

        % Weighting factors
        f2 = d2/(d2+d1+dD+dbp);
        f1 = d1/(d2+d1+dD+dbp);
        fD = dD/(d2+d1+dD+dbp);
        fbp = dbp/(d2+d1+dD+dbp);
        
        % Effective Thermal Properties
        tceffb = ((f2/tc2)+(f1/tc1)+(fD/tcD)+(fbp/tcbp) )^-1;
        vhceffb =  (f2*vhc2)+(f1*vhc1)+(fD*vhcD)+(fbp*vhcbp) ;
        tdeffb = tceffb/vhceffb;
        
        hpf = 0.5;   % Heat partition factor (–)
      
        xspace = linspace(-5*b,2*b,xsteps); % Position vector for distance along travel direction  
        k = (v/(2*tdeffb));                 % Convienent grouping of parameters
        pre = (2*sqrt(2)*freq*cof*amp*hpf*load)/(pi*tceffb*sw); % Convienent grouping of parameters
        
        % Defining Temperature Functions
        f = @(X,x) (pre/b)*real(sqrt(1 - (X/b).^2)).*exp(-k.*(x-X)).*besselk(0,k.*sqrt((x-X).^2));
        g = @(x) integral(@(X) f(X,x), -b,b); 
            for xi = 1:xsteps
                temp(xi) = tamb + g(xspace(xi));
            end
           
            time = ((xspace)/v) - ((xspace(1))/v);              %   Populates time vector
            dif =  DL.* exp(-QL./((8.617e-5)*temp));            %   Lattice interdiffusivity
            difgb = DGB.* exp(-QGB./((8.617e-5)*temp));         %   Grain boundary interdiffusivity
            difeff = dif + ((3*0.6)/15)*difgb;                  %   Effective interdiffusivity
            dif = difeff;
            
            difd =  2*max(sqrt(cumtrapz(time, dif)));           %   Estimate for diffusion distance
            leadedgetime = (b/v);                               %   Half sonotrode contact time
            result = find(time < leadedgetime);
            dif(result) = 0;
            
            % Crank-Nicolson Solver 
            dx=difd/100;   
            nx=500;     
            nt=500;        
            zv = dx*linspace(-nx/2,nx/2,nx);
            dt = time(55)-time(54);
            tv = dt*linspace(1,nt,nt);
            xv = ((tv*v))/1e-3;
              
                UUU = zeros(nt,nx);
                Uo(1:nx/2)=1; Uo(nx/2+1:nx)=0; 
                Un(nx)=0; 
                Un(1)=1;
                Un(nx/2)=0.5; 
                UUU(1,:)=Uo;
                for k=2:nt
                    alfa = dif(k-1);
                    b1 = alfa/(2*dx^2);
                    c = b1; 
                    a = 1/dt+b1+c; 

                    for ii=1:nx-2
                        if ii==1
                            d(ii)=c*Uo(ii)+(1/dt-b1-c)*Uo(ii+1)+b1*Uo(ii+2)+c*Un(1);
                        elseif ii==nx-2
                            d(ii)=c*Uo(ii)+(1/dt-b1-c)*Uo(ii+1)+b1*Uo(ii+2)+b1*Un(nx);
                        else
                            d(ii)=c*Uo(ii)+(1/dt-b1-c)*Uo(ii+1)+b1*Uo(ii+2);
                        end
                    end 
                 bb=b1*ones(nx-3,1);
                 cc=bb;
                 aa=a*ones(nx-2,1);
                 AA=diag(aa)+ diag(-bb,1)+ diag(-cc,-1);
                 UU=AA\d';
                 Un=[Un(1),UU',Un(nx)];
                 UUU(k,:)=Un;
                 Uo=Un;
                end

                UUU = rot90(UUU);
                UUU = fliplr(UUU);
                [c, index] = min(abs(UUU(:,1)-0.01));     
                
                % Populates diffusion distance matrix
                dmat(int16(load*lstep/llim),int16(v*vstep/vlim)) = 2*abs(zv(index));
               
    end
    progress = load*lstep/llim
end
close all


dcon = [0.01,0.05,0.1,0.2, 0.5, 1];    % Chose diffusion distance contours (um)

[c,w] = contour(vv,lv,dmat,dcon);  
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
    
mat = zeros(1000,2*length(dcon)+4);

for k = 1:length(dcon)
    x = linspace(s(k).xdata(1), s(k).xdata(end), 1000);
    y = interp1(s(k).xdata, s(k).ydata, x);
    mat(:, 2*k -1) = x;
    mat(:,2*k) = y;
end

close all

% Plot diffusion distance contours
figure
plt = plot(mat(:,1),mat(:,2), mat(:,3),mat(:,4),mat(:,5),...
    mat(:,6),mat(:,7),mat(:,8),mat(:,9),...
    mat(:,10),mat(:,11),mat(:,12));
set(plt,'LineWidth',2);    
set(gca,'linewidth',1.5)
set(gca, 'FontSize', 16)
xlabel('sonotrode velocity', 'FontSize',20)
ylabel('weld pressure (MPa)', 'FontSize',20)





