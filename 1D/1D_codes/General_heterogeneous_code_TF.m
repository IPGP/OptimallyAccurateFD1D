% General code for 1D homogeneus & isotrope elastic case without internal
% boundaries & without PML

%

%  The original code is written by NF, OO, CC
%   

%
%             Here're Thijs Franken and Guy Mabylaht working (Summer 2016)




close all ;
clear all ;
%clc ;

%-------------------------------------------------------------------------
%                        MODELISATION SETTINGS
%-------------------------------------------------------------------------

SNAPSHOT = false ; % Makes snapshots of the modelisation
Snapshot_Time = 340:20:600 ; % Parameters of snapshots temporisation

MAKE_MOVIE = false ; % Saves a video of the modelisation
tag = 'General 1D' ;

SAVE_SEISMOGRAMS = false ; % Saves a final steps snapshot of the seismograms
seis_tag = ['General_1D_Code'] ;

choice = menu('Please choose a method','CONV 2','CONV 4', 'OPT 2', 'OPT 4') ;

if choice == 1 ;
    CONV2 = true ;
    OPT2 = false ;
    CONV4 = false ;
    OPT4 = false ;
elseif choice == 2 ;
    CONV2 = false ;
    OPT2 = false ;
    CONV4 = true ;
    OPT4 = false ;
elseif choice == 3 ;
    CONV2 = false ;
    OPT2 = true ;
    CONV4 = false ;
    OPT4 = false ;
else
    CONV2 = false ;
    OPT2 = false ;
    CONV4 = false ;
    OPT4 = true ;
end

% CONV2 = true ;
% CONV4 = false ;
% OPT2 = false ;
% OPT4 = false ;

prompt = {'Amount of nodes','Distance of propagation [m]',...
          'P-Waves velocity model','Density model',...
          'Size of the time steps [ms]','Frequency of the source [Hz]',...
          'Source positions [m]','Amount of receivers'};
settitle = '1D-Wave modelisation settings';
num_lines = 1 ;
usual_values = {'200','100','9999','9999','0.1','100','50','2'} ;
answer = inputdlg(prompt,settitle,num_lines,usual_values) ;

set = zeros(8,1);
for i = 1:8
    set(i,1) = str2num(answer{i});
end

if set(7) >= set(2)
    fprintf('Please choose a source positions within propagation domains. \n')
    return ;
end

input('Start time loop?')
tic ;

%-------------------------------------------------------------------------
%                          DEF OF THE 1D CASE
%-------------------------------------------------------------------------


% Grid's definition (space, properties & time)
%-------------------------------------------------------------------------

%NX = input('Amount of nodes : ') ;
%NX = 200.d0 ; % Nb of nodes
NX = set(1) ;
mg = NX / 100.d0 ;

%XMAX = input('Distance of propagation [m] : ') ;
%XMAX = 100.d0 ; % Maximal X value [m]
XMAX = set(2) ;
XMIN = 0.d0 ; % Minimal X value [m]

DELTAX = (XMAX-XMIN) / NX ; % Space between nodes [m]
dx2 = DELTAX^2.d0 ;

xvec = [0:NX] * DELTAX ; % Vector containing all of the X steps

%cp = input('P-Waves velocity [km/s] : ') ;
%cp = 1800.d0 ; % Velocity of P-Waves in the medium [km/s]
cpmodel = set(3) ;
%cs = input('S-Waves velocity [km/s] : ') ;



%rho = 2400.d0 ; % Density of the medium [-]
rhomodel = set(4) ;




%definition of the variables for holding the medium parameters
rho=zeros(NX+1,1);  %density
mu=zeros(NX+1,1); % stiffness 
cp=zeros(NX+1,1); % p wave velocity
cs=zeros(NX+1,1); % s wave velocity


% initialization of density and P-wave velocity for each node
rho(:,1)=(linspace(2400,2500,NX+1));
cp(:,1)=(linspace(1500,2000,NX+1)); %[km/s]
cs(:,1)=cp(:,1)/1.732d0; %[km/s]
 
%computation of lame parameter mu
mu(:,1)=rho(:,1).*cs(:,1).*cs(:,1);


%fprintf('SPACE:\n NX = %d \n XMAX = %.2f \n XMIN = %.2f \n dX = %.2f\n', NX, XMAX, XMIN, DELTAX) ;
%fprintf('PROPERTIES:\n cp = %.0f \n cs = %.0f \n rho = %.0f \n lambda = %d\n mu = %d\n', cp, cs, rho, lambda, mu) ;

%fprintf('PROPERTIES:\n cp = %.0f \n cs = %.0f \n rho = %.0f \n mu = %d\n', cp, cs, rho, mu) ;


%----------

%DELTAT = input('Size of the time steps [s] : ') ;
%DELTAT = 1.d-4 ; % Size of time steps [s]
DELTAT = set(5) ;
DELTAT = DELTAT*(10^-3);
DELTAT = (DELTAT) / mg ;
dt2 = DELTAT^2.d0 ;
time = 0.08d0 ;

NSTEP = round(time/DELTAT) ; % Nb of time steps [s]

time_vec = [1:NSTEP]' * DELTAT ; % Vector containing all of the time steps


IT_DISPLAY = 50 ; % Nb of iterations between each display

fprintf('TIME:\n IT_DISP = %d \n NSTEP = %d \n dT = %d\n', IT_DISPLAY, NSTEP, DELTAT) ;


% Def of other variables & matrices
%--------------------------------------------------------------------------

ZERO = 0.d0 ; % Def of zero

e = 0.00001d0 ; % Def of an infinitesimal nb (~zero)
E = 1.d+30 ;  % Def of a huge nb (~infinite)

STABILITY_THRESHOLD = 1.d+25 ; % Def of a Vmax. If Vp&Vs>Vmax : Unstable code


% Def of the source
%--------------------------------------------------------------------------

%F0 = input('Frequency of the source [Hz] : ');
%F0 = 100.d0 ; % Frequency [Hz]
F0 = set(6) ;
T0 = 1.20d0 / F0 ; % Time of source triggering [s]

Ampl = 1.d9 ; % Amplitude of the signal

Source_Pos = set(7) ; % Source position [m]
mindist = XMAX ;

% Way to determinate the closest-to-source node
for i = 2:NX
    
    dist = Source_Pos - DELTAX*i ;
    
    if abs(dist) < mindist
        mindist = dist ;
        ISOURCE = i ;
    end
    
end

clearvars mindist dist ;

fprintf('Source position: %.2f m i=%d\n', Source_Pos, ISOURCE) ;


% Def of the receivers
%--------------------------------------------------------------------------

%NREC = input('Choose the number of receivers : ');
%NREC = 21 ; % Nb of receivers
NREC = set(8) ;

xdeb = 0.1d0 * XMAX ; % Position of the first receiver [m]
xfin = 0.9d0 * XMAX ; % Position of the last receiver [m]

xrec = zeros(NREC,1) ; % Def vector containing the receivers positions
ix_rec = zeros(NREC,1) ; % Def vector containing receivers closest-nodes positions

xspacerec = (xfin-xdeb) / double(NREC-1) ; % Receivers spacing

% Way to "place" receivers along the grid
for irec = 1:NREC
    xrec(irec) = xdeb + double(irec-1) * xspacerec ;
end

% Way to find the closest node of each receiver
% (Same way as to find the closest-to-source node)
for irec = 1:NREC
    
   dist = E ;
   
   for i = 1:NX
       
       distval = abs(DELTAX*double(i-1) - xrec(irec)) ;
       
       if (distval < dist)
           dist = distval ;
           ix_rec(irec) = i ;
       end
       
   end
   
   fprintf('Receiver %d x= %.2f\n',irec,xrec(irec))
   fprintf('Closest node found at %.2f m from i = %d\n\n',dist,ix_rec(irec));
   
end

%-------------------------------------------------------------------------
%                 CHECK OF THE MODELISATIONS CONSISTENCY
%-------------------------------------------------------------------------
  
Courant_number = max(cp)* DELTAT * DELTAX ;

fprintf('Courant Number = %.4f\n', Courant_number) ;

if Courant_number > 1.d0
    disp('Time step is too large, the simulation will be unstable. \n Please choose another time step.') ;
    return ;
end


%--------------------------------------------------------------------------
%                    MODELISATIONS OBJECTS PARAMETERS
%--------------------------------------------------------------------------


% Initiation of a videos object
if MAKE_MOVIE
    movie_name_vy = ['1D_u_' num2str(NX) '_' num2str(DELTAX) '_' num2str(F0) '.avi'] ;
	vidObj_vy = VideoWriter(movie_name_vy) ;
	open(vidObj_vy) ;
end

  

%--------------------------------------------------------------------------
%                          MODELISATIONS CODE
%--------------------------------------------------------------------------

a = pi*pi * F0*F0 ;
umaxv = Ampl * exp(-a*T0^2) / 2.d0 * DELTAT ; % Def of a display maximum

u = zeros(NX+1,1) ;
unm1 = zeros(NX+1,1) ;
unm2 = zeros(NX+1,1) ;
% 
% u1 = zeros(NX+1,1) ;
% u1nm1 = zeros(NX+1,1) ;
% u1nm2 = zeros(NX+1,1) ;

UI2 = zeros(3,3) ;
UI4 = zeros(5,3) ;
g = zeros(NX+1,1) ;
utotal = zeros(NX+1,1) ;

A_con2 = zeros(3,3,NX-1) ;
K_con2 = zeros(3,3,NX-1) ;
A_opt2 = zeros(3,3,NX-1) ;
K_opt2 = zeros(3,3,NX-1) ;
dA2 = zeros(3,3,NX-1) ;
dK2 = zeros(3,3,NX-1) ;





 for i=2:NX
    
    A_con2(1:3,1:3,i) = (1/dt2)*[0, rho(i), 0;0,-2*rho(i),0;0,rho(i),0];
    
    K_con2(1:3,1:3,i) = (1/2/dx2)*[0, 0 ,0; (mu(i-1,1)+mu(i,1)), -(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)), (mu(i,1)+mu(i+1,1)); 0 0 0];
    
    A_opt2(1:3,1:3,i) = (1/(12*dt2))*[rho(i), 10*rho(i), rho(i); -2*rho(i), -20*rho(i),-2*rho(i); rho(i), 10*rho(i), rho(i)];
    
    K_opt2(1:3,1:3,i)= (1/(24*dx2))*[(mu(i-1,1)+mu(i,1)), -(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)), (mu(i,1)+mu(i+1,1));...
        10*(mu(i-1,1)+mu(i,1)), -10*(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)), 10*(mu(i,1)+mu(i+1,1));...
        (mu(i-1,1)+mu(i,1)), -(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)), (mu(i,1)+mu(i+1,1))]; 
    
    dA2(1:3,1:3,i)= A_opt2(1:3,1:3,i)-A_con2(1:3,1:3,i);
   
    dK2(1:3,1:3,i)= K_opt2(1:3,1:3,i)-K_con2(1:3,1:3,i);
  
  end  

% 
% 
% 
% for i = 2:NX
%   A_con2(1:3,1:3,i-1) = (rho/dt2)*[0 1 0; 0 -2 0; 0 1 0] ;
%   K_con2(1:3,1:3,i-1) = (mu/dx2)*[0 0 0; 1 -2 1; 0 0 0] ;
%   A_opt2(1:3,1:3,i-1) = (rho/(12*dt2))*[1 10 1; -2 -20 -2; 1 10 1] ;
%   K_opt2(1:3,1:3,i-1) = (mu/(12*dx2))*[1 -2 1; 10 -20 10; 1 -2 1] ;
%   dA2(1:3,1:3,i-1) = (rho/(12*dt2))*[1 10 1; -2 -20 -2; 1 10 1]-(rho/dt2)*[0 1 0; 0 -2 0; 0 1 0] ;
%   dK2(1:3,1:3,i-1) = (mu/(12*dx2))*[1 -2 1; 10 -20 10; 1 -2 1]-(mu/dx2)*[0 0 0; 1 -2 1; 0 0 0] ;
% end

A_con4 = zeros(3,5,NX-1) ;
K_con4 = zeros(3,5,NX-1) ;
A_opt4 = zeros(3,5,NX-1) ;
K_opt4 = zeros(3,5,NX-1) ;
dA4 = zeros(3,5,NX-1) ;
dK4 = zeros(3,5,NX-1) ;



 for i=3:NX-1
    
     A_con4(1:3,1:5,i) = (1/dt2)*[0, 0, rho(i), 0, 0; 0, 0, -2*rho(i), 0, 0; 0, 0, rho(i), 0, 0];
     
     K_con4 (1:3,1:5,i)= (1/(24*dt2))*[0 ,0,0,0,0;...
         -(mu(i-2,1)+mu(i,1)), 16*(mu(i-1,1)+mu(i,1)), -16*(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)),16*(mu(i,1)+mu(i+1,1)), -(mu(i,1)+mu(i+2,1))+(mu(i-2,1)+2*mu(i,1)+mu(i+2,1));...
         0,0,0,0,0];
     
     A_opt4(1:3,1:5,i)= (rho(i)/(90*dt2))*[-1 4 84 4 -1 ; 2 -8 -168 -8 2 ; -1 4 84 4 -1];
     
    
     K_opt4(1:3,1:5,i)=(1/(144*dt2))*[-(mu(i-2,1)+mu(i,1)), 16*(mu(i-1,1)+mu(i,1)), -16*(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)),16*(mu(i,1)+mu(i+1,1)), -(mu(i,1)+mu(i+2,1))+(mu(i-2,1)+2*mu(i,1)+mu(i+2,1));...
         -10*(mu(i-2,1)+mu(i,1)), 16*10*(mu(i-1,1)+mu(i,1)), -16*10*(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)),16*10*(mu(i,1)+mu(i+1,1)), -10*((mu(i,1)+mu(i+2,1))+(mu(i-2,1)+2*mu(i,1)+mu(i+2,1)));...
         -(mu(i-2,1)+mu(i,1)), 16*(mu(i-1,1)+mu(i,1)), -16*(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)),16*(mu(i,1)+mu(i+1,1)), -(mu(i,1)+mu(i+2,1))+(mu(i-2,1)+2*mu(i,1)+mu(i+2,1))];
          
    
    dA4(1:3,1:5,i)= A_opt4(1:3,1:5,i)-A_con4(1:3,1:5,i);
   
    dK4(1:3,1:5,i)= K_opt4(1:3,1:5,i)-K_con4(1:3,1:5,i);
  
  end  




% 
% for i = 3:NX
%    A_con4(1:3,1:5,i-2) = (rho/dt2)*[0 0 1 0 0; 0 0 -2 0 0; 0 0 1 0 0] ;
%    K_con4(1:3,1:5,i-2) = (mu/(12*dx2))*[0 0 0 0 0 ; -1 16 -30 16 -1 ; 0 0 0 0 0] ;
%    A_opt4(1:3,1:5,i-2) = (rho/(90*dt2))*[-1 4 84 4 -1 ; 2 -8 -168 -8 2 ; -1 4 84 4 -1] ;
%    K_opt4(1:3,1:5,i-2) = (mu/(144*dx2))*[-1 16 -30 16 -1; -10 160 -300 160 -10;-1 16 -30 16 -1] ;
%    dA4(1:3,1:5,i-2) = (rho/(90*dt2))*[-1 4 -6 4 1; 2 -8 12 -8 2 ;-1 4 -6 4 1] ;
%    dK4(1:3,1:5,i-2) = (mu/(144*dx2))*[-1 16 -30 16 -1; 2 -32 60 -32 2;-1 16 -30 16 -1] ;
% end
% 

%--------------------------------------------------------------------------


for it = 1:NSTEP
    
    t = double(it-1) * DELTAT ;
    source_term = Ampl * exp(-a*(t-T0)^2)/ rho(ISOURCE) ;
      
% CONV2
%--------------------------------------------------------------------------
    if CONV2
           
       PAUSE_ON = false ; % Adds a pause between each snapshot
       Pause_Time = 0.5 ; % Duration of that pause [s]
          
       u(:) = 0.d0 ;
          
       for i = 2:NX-2
           value_mu_du_dxx = unm1(i-1)*K_con2(2,1,i)+unm1(i)*K_con2(2,2,i)+unm1(i+1)*K_con2(2,3,i); 
           u(i) = (value_mu_du_dxx-unm1(i)*A_con2(2,2,i)-unm2(i)*A_con2(3,2,i))/A_con2(1,2,i);
       end
      
       u(ISOURCE) = source_term * DELTAT ;

       % Time steps update
       unm2 = unm1 ;
       unm1 = u ;
          
          
% OPT2
%--------------------------------------------------------------------------
    elseif OPT2
           
       PAUSE_ON = false ;
           
       u(:) = 0.d0 ;
       g(:)= 0.d0 ;
    
       
       % Predictor
       
       for i= 2:NX
           value_mu_du_dxx = unm1(i-1)*K_con2(2,1,i)+unm1(i)*K_con2(2,2,i)+unm1(i+1)*K_con2(2,3,i);
           u(i) = (value_mu_du_dxx-unm1(i)*A_con2(2,2,i)-unm2(i)*A_con2(3,2,i))/A_con2(1,2,i);
       end
    
          
       g(:) = 0.d0 ;
          
       
       % Corrector
       
       for i = 2:NX
           UI2(1:3,1:3,i-1) = [u(i-1) u(i) u(i+1); ...
                              unm1(i-1) unm1(i) unm1(i+1);...
                              unm2(i-1) unm2(i) unm2(i+1)] ;
           g(i) = -sum(sum((dA2(1:3,1:3,i)-dK2(1:3,1:3,i)).*UI2(1:3,1:3,i-1)))/A_con2(1,2,i) ;
       end
          
       for i = 2:NX
           u(i) = u(i)+(g(i)/A_con2(1,2,i)); %update: c = c0 + kdc
       end
          
       u(ISOURCE) = source_term * DELTAT ;
          
       % Time steps update
       unm2 = unm1 ;
       unm1 = u ;
          
    
% CONV4
%--------------------------------------------------------------------------
    elseif CONV4
           
        PAUSE_ON = false ; % Adds a pause between each snapshot
        Pause_Time = 0.1 ; % Duration of that pause [s]
                      
        u(:) = ZERO ;
           
        for i = 3:(NX-1)
            value_mu_du_dxx = unm1(i-2)*K_con4(2,1,i)+unm1(i-1)*K_con4(2,2,i)+unm1(i)*K_con4(2,3,i)+...
                              unm1(i+1)*K_con4(2,4,i)+unm1(i+2)*K_con4(2,5,i);
                          
            u(i) = (value_mu_du_dxx-unm1(i)*A_con4(2,3,i)-unm2(i)*A_con4(3,3,i))/A_con4(1,3,i);
        end
      
        i = ISOURCE ;
        u(i) = source_term * DELTAT ;
          
        % Time steps update
        unm2 = unm1 ;
        unm1 = u ;
           
           
% OPT4
%--------------------------------------------------------------------------           
    elseif OPT4
           
        PAUSE_ON = false ;
          
        u(:) = 0.d0 ;
        g(:)= 0.d0 ;
        
   
        utotal(:) = 0.d0;
        UI4(:,:) = 0.d0;
           
        
        % Predictor
           
        for i = 3:(NX-1)
            value_mu_du_dxx = unm1(i-2)*K_con4(2,1,i)+unm1(i-1)*K_con4(2,2,i)+unm1(i)*K_con4(2,3,i)+...
                              unm1(i+1)*K_con4(2,4,i)+unm1(i+2)*K_con4(2,5,i);
            u(i) = (value_mu_du_dxx-unm1(i)*A_con4(2,3,i)-unm2(i)*A_con4(3,3,i))/A_con4(1,3,i);
        end
           
        g(:)=0.d0;
           
        
        % Corrector
        
        for i = 3:(NX-1)
            UI4(1:3,1:5,i-2) = [u(i-2) u(i-1) u(i) u(i+1) u(i+2) ; ...
                   unm1(i-2) unm1(i-1) unm1(i) unm1(i+1) unm1(i+2); ...
                   unm2(i-2) unm2(i-1) unm2(i) unm2(i+1) unm2(i+2)];
            g(i) = -sum(sum((dA4(1:3,1:5,i)-dK4(1:3,1:5,i)).*UI4(1:3,1:5,i-2)))/A_con4(1,3,i);
        end
           
        for i=3:(NX-1)
            u(i) = u(i)+(g(i)/A_con4(1,3,i)); %update: c = c0 + kdc
        end
           
        u(ISOURCE) = (source_term * DELTAT) ;  
           
        % Time steps update
        unm2 = unm1;
        unm1 = u;
        
           
    end
       
       
%--------------------------------------------------------------------------

       % Dirichlet Boudaries Conditions
%        u(1) = ZERO ;
%        u(NX+1) = ZERO ;
           
  
% DISPLAY PART
%--------------------------------------------------------------------------

    if SAVE_SEISMOGRAMS
        for irec = 1:NREC
            seisux(it,irec) = u(ix_rec(irec));
        end 
    end
        
    if(mod(it,IT_DISPLAY) == 0)
        
        
        if max(u) > umaxv
            umaxv = max(u) ;
        end
        
        plot(xvec,u,'linewidth',2, 'Color','r') ;
       
        grid on ;
        axis([min(xvec) max(xvec) -umaxv umaxv]) ;
        xlabel('x, m','fontSize',14) ;
        ylabel('u','fontSize',14) ;
        titlestring = ['TIME STEP = ',num2str(it), ' TIME = ',num2str(t), ' s'] ;
        title(titlestring ,'fontsize',14) ;
        drawnow ;        
        
%         if u(1) || u(NX) > 1.d-3
%             toc ;
%             return ;
%             fprintf('End.')
%         end
        
        if PAUSE_ON
            pause(Pause_Time);
        end
        
        if MAKE_MOVIE
            F_y = getframe(gcf) ;  %-  capture figure or use gcf to get current figure. Or get current
            writeVideo(vidObj_vy,F_y);  %- add frame to the movie
            %fprintf('Frame for %s captured\n',movie_name_vy);
        end
            
        if SNAPSHOT
            if  nnz(snapshot_time==it)>0
                snapshat = getframe(gcf);
                imgg = frame2im(snapshat);
                scrsht_name=['im' num2str(it) '.png'];
                imwrite(imgg,scrsht_name);
                fprintf('Screenshot %s saved to %s\n', scrsht_name, pwd);
                clearvars scrsht_name imgg snapshat
            end  
        end
        
    end
end



%--------------------------------------------------------------------------
%                         MODELISATIONS RECORD
%--------------------------------------------------------------------------

if SAVE_SEISMOGRAMS
    for i = 1:NREC
        filename = [seis_tag 'x' num2str(xrec(i),'%.2f') 'seis_' num2str(i) '.txt'] ;
        dlmwrite(filename, [time_vec, seisux(:,i)]) ;
        fprintf('Seismogram for rec at %.2f saved as %s to %s\n', xrec(i), filename, pwd) ;
    end
end

if MAKE_MOVIE
    close(vidObj_vy) ;
    fprintf('Video %s saved in %s\n',movie_name_vy, pwd) ;
end

toc ;

disp('End.') ;