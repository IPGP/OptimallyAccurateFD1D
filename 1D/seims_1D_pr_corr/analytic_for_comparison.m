% SOLUTION Analytique 

% 1D case disp

close all;
clear all;
clc;
 mg=2;
 NX =100*mg;  %X
 
 time=0.08d0;
 
 % time step in seconds
 DELTAT = 1.d-4;	%[sec] 
 DELTAT = DELTAT/mg;
 % total number of time steps
 %  NSTEP = 2000;
 NSTEP = round(time/DELTAT);
 time_vec = [1:NSTEP]'*DELTAT;
 
XMAX=100;
XMIN=0.d0;
XFE=XMAX/3.d0;

IT_DISPLAY = 50;
   
% NSTEP=1000;
% DELTAT = 1.d-4;	%[sec]

%Pause a little bit each iteration
   PAUSE_ON=true;
   pause_time=0.1; %[sec]
   
FE_BOUNDARY=false; %just adds Cijkl step
MAKE_MOVIE=false;
DATA_TO_BINARY_FILE=false;
SNAPSHOT=false;
snapshot_time=340:20:600; %on what steps
tag='pure';

SAVE_SEISMOGRAMS=false;
% seis_tag='mz200';
seis_tag=['pure1D' num2str(NX)];

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

DELTAX=(XMAX-XMIN)/NX;

% parameters for the source

f0 = 100.d0;
t0 = 1.20d0 / f0;
factor = 1.d9;
a = pi*pi*f0*f0;

cp_above_eb=1800.d0;
 cp_below_eb=1800.d0;
 rho_above_eb=2400.d0;
 rho_below_eb=2400.d0;
 
 
 
 
order=1;
C=zeros(NX,order,order,order,order);
densitya = rho_above_eb;
cpa = cp_above_eb;      %[km/s]
csa = cpa / 1.732d0;	%[km/s]
lambdaa =densitya*(cpa*cpa - 2.d0*csa*csa);
mua = densitya*csa*csa;
densityb = rho_below_eb;
cpb = cp_below_eb;	%[km/s]
csb = cpb / 1.732d0;	%[km/s]
lambdab =densityb*(cpb*cpb - 2.d0*csb*csb);
mub = densityb*csb*csb;
% 
% 
% 
% 
% 
% % omega - k formulation
% 
% Nomega=2^11;
% Nk=2^11;
% u_k=zeros(2*Nk,2*Nomega);
% u=zeros(2*Nk,2*Nomega);
% for ik=1:Nk;
%     u_freq=zeros(2*Nomega,1);
%     for iomega=1:Nomega;
%         omega=2*pi/time*iomega;
%         k=2*pi/XMAX;
%         u_freq(iomega)=1/(mua*k*k-densitya*omega*omega);
%     end
%     u_k(ik,:)=ifft(u_freq);
% end
% 
% for iomega=1:2*Nomega;
%     u(:,iomega)=ifft((u_k(:,iomega)));
% end
% 
% 
% 
% break;
% 
% 
% 
% 
% 
% 
% 
% % omega - x formulation
% 
% Nomega=2^11; 
% u_freq=zeros(2*Nomega,1);
% 
% for i=1:NX+1
%     x_trial=(i-1)*DELTAX;    
%     for iomega=1:Nomega;
%         omega=2*pi/time*iomega;
%         u_freq(iomega)=-1i*csa/(4*pi*mua)/omega*exp(-1i*omega/csa*x_trial);
%     end
%     u(i,:)=ifft(u_freq);
% end

    


for j=1:NSTEP
    
    t = double(j-1)*DELTAT;
    %source_term = factor*exp(-a*(t-t0)^2);
    source_term=0.;
    if t == 0.
        source_term=1.;
    end
        
    ricker = source_term*DELTAT^2/densitya;
    
    for i = 2:NX
        
        step = i*DELTAX;
        c = csa*j*DELTAT^2/densitya;
        
        if step>(-c) && step>c
            u(i)= conv(ricker,pi/(densitya*csa));
        else
            u(i)=conv(ricker,0);
        end
    end
end

u(1)=0.d0;
u(NX+1)=0.d0;
