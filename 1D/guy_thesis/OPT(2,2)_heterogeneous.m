% Optimally accurate finite difference operators, second order in time and
% fourth order in space. 


 mg=2;   % factor for defining the number of Grid points
 NX =100*mg;  % number of grid points in the model

 
 % time step in seconds
  time=0.08d0;
 DELTAT = 1.d-4;	%[sec] 
 DELTAT = DELTAT/mg;
 
 % total number of time steps
 NSTEP = round(time/DELTAT);  % computes the # of times steps in 0.08 second
 time_vec = (1:NSTEP)'*DELTAT; % stores the time steps into an array
 
%definition of the model size 
XMAX=100; % length of the model 
XMIN=0.d0;
XFE=XMAX/3.d0;
DELTAX=(XMAX-XMIN)/NX;  % spacing between the nodes
xvec=(0:NX)*DELTAX; % definition of the Grid points with their spacing

IT_DISPLAY = 10;
   

%Pause a little bit each iteration
   PAUSE_ON=true;
   pause_time=0.01; %[sec]
   
FE_BOUNDARY=false; %just adds Cijkl step
MAKE_MOVIE=false;
DATA_TO_BINARY_FILE=false;
SNAPSHOT=false;
snapshot_time=340:20:600; %on what steps
tag='pure';

SAVE_SEISMOGRAMS=false;
% seis_tag='mz200';
seis_tag=['pure1D' num2str(NX)];


% parameters for the source
f0 = 100.d0;
t0 = 1.20d0 / f0;
factor = 1.d9;

%Define i of source
xsource=0.01d0*XMAX;    %arbitrary source position in [m]
mindist=XMAX;           
for i=3:NX
    dist=xsource-DELTAX*i;
    if abs(dist)<mindist
        mindist=dist;
        ISOURCE=i;
    end
end
clearvars mindist dist;

%displays information about the Space Grid and Time Grid
fprintf('TIME:\n IT_DISP=%d NSTEP=%d\n', IT_DISPLAY, NSTEP);
fprintf('SPACE:\n NX=%d XMAX=%.2f XMIN=%.2f dX=%.2f\n', NX, XMAX, XMIN, DELTAX);
fprintf(' FE Boundary at XFE=%.2f\n', XFE);

% zero
 ZERO = 0.d0;
 eps=0.00001d0;

% large value for maximum
 HUGEVAL = 1.d+30;
 

%definition of the variables for holding the medium parameters
rho=zeros(NX+1,1);  %density
mu=zeros(NX+1,1); % stiffness 
cp=zeros(NX+1,1); % p wave velocity
cs=zeros(NX+1,1); % s wave velocity
 

% Acquisition parameters
NREC=21;  % number of receivers
xdeb=0.2d0*XMAX; % left end point of the receiver line
xfin=0.6d0*XMAX; % right end point of the receiver line
xrec=zeros(NREC,1); % receiver line definition
ix_rec=zeros(NREC,1);
xspacerec = (xfin-xdeb) / double(NREC-1); % (xfin - xdeb) difference between left and right end points (length of the receiver line), xspacerec is the receiver spacing

% initialization of the receiver line (receivers are placed at their
% respective positions (20, 20+xspacerec, 20+2*(xspacerec), etc)
for irec=1:NREC
    xrec(irec) = xdeb + double(irec-1)*xspacerec;
end


% find closest grid point for each receiver
for irec=1:NREC
   dist = HUGEVAL;
    for i = 1:NX
      distval = abs(DELTAX*double(i-1) - xrec(irec));
      if(distval < dist)
        dist = distval;
        ix_rec(irec) = i;
      end
   end
   fprintf('receiver %d x= %.2f\n',irec,xrec(irec));
   fprintf('closest grid point found at distance %.2f in i = %d\n\n',dist,ix_rec(irec));
end


% initialization of density and P-wave velocity for each node
rho(:,1)=(linspace(2400,2500,NX+1));
cp(:,1)=(linspace(1500,2000,NX+1)); %[km/s]
cs(:,1)=cp(:,1)/1.732d0; %[km/s]
 
%computation of lame parameter mu
mu(:,1)=rho(:,1).*cs(:,1).*cs(:,1);

% definition Cijkl matrix
order=1;
C_ijkl=zeros(NX+1,1);
C_ijkl(:,1)=mu(:,1);


Courant_number1 = 2400 * DELTAT * DELTAX;
  fprintf('Courant number + = %.4f\n',Courant_number1);
  
  Courant_number2 = 2480 * DELTAT * DELTAX;
  fprintf('Courant number - = %.4f\n\n',Courant_number2);

  if(Courant_number1 > 1.d0) || (Courant_number2 > 1.d0) 
      disp('time step is too large, simulation will be unstable');
  end
  
  
  %source parameters
   a = pi*pi*f0*f0;
  umaxv = factor*exp(-a*t0^2)/2.d0*DELTAT;
  
  % delta x and t to the second
  dx2 = DELTAX^2.d0;
  dt2 = DELTAT^2.d0;



  
   % definition of the displacement variables 
u=zeros(NX+1,1);
u1=zeros(NX+1,1);
unm1=zeros(NX+1,1);
unm2=zeros(NX+1,1);
u1nm1=zeros(NX+1,1);
u1nm2=zeros(NX+1,1);
utotal=zeros(NX+1,1);
  UI=zeros(3,3);
  fs=zeros(3,3);
  f1s=zeros(3,3);
  fs_corr=zeros(3,3);
  g=zeros(NX+1,1);
  g0=zeros(NX+1,1);
  g1=zeros(NX+1,1);
  
  %definition of conventional operators in matricial form
  A_con=zeros(3,3,NX-1);
  K_con=zeros(3,3,NX-1);
  
  %definition of optimally accurate operators in matricial form (will not be used in the
  %scheme because they yield an implicit scheme)
  A_opt=zeros(3,3,NX-1);
  K_opt=zeros(3,3,NX-1);
  
  %error or residuals between conventional and optimally accurate operators
  %in matricial form
  dA=zeros(3,3,NX-1);
  dK=zeros(3,3,NX-1);
  
  input('Start time loop?')
  
  
  %computation of Optimally accurate and conventional matricial
  %operators as well as of their corresponding residuals
  
  for i=2:NX-1
    
    A_con(1:3,1:3,i) = (1/dt2)*[0, rho(i,1), 0;0,-2*rho(i),0;0,rho(i),0];
    
    K_con(1:3,1:3,i) = (1/dx2)*[0, 0 ,0; (mu(i-1,1)+mu(i,1)), -(mu(i-1,1)+2*mu(i,1)+mu(i,1)+mu(i+1,1)), (mu(i,1)+mu(i+1,1)); 0 0 0];
    
    A_opt(1:3,1:3,i) = (1/(12*dt2))*[rho(i), 10*rho(i), rho(i); -2*rho(i), -20*rho(i),-2*rho(i); rho(i), 10*rho(i), rho(i)];
    
    K_opt(1:3,1:3,i)= (1/(24*dx2))*[(mu(i-1,1)+mu(i,1)), -(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)), (mu(i,1)+mu(i+1,1));...
        10*(mu(i-1,1)+mu(i,1)), -10*(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)), 10*(mu(i,1)+mu(i+1,1));...
        (mu(i-1,1)+mu(i,1)), -(mu(i-1,1)+2*mu(i,1)+mu(i+1,1)), (mu(i,1)+mu(i+1,1))]; 
    
    dA(1:3,1:3,i)= A_opt(1:3,1:3,i)-A_con(1:3,1:3,i);
   
    dK(1:3,1:3,i)= K_opt(1:3,1:3,i)-K_con(1:3,1:3,i);
  
  end  
  
  
  
  
  % Predictor-corrector scheme
  
  for it = 1:NSTEP
        
    u(:)=0.d0;
    u1(:)=0.d0;
    UI(:,:)=0.d0;

   
    % Ricker source time function (second derivative of a Gaussian)
     t = double(it-1)*DELTAT;
    source_term = factor * exp(-a*(t-t0)^2);
    
    % Predictor
    g(:)=0.d0;
    g(ISOURCE)=source_term*DELTAT^2/rho(ISOURCE); 
        
    
        %wavefield extrapolation with conventional operators
    for i=2:NX
        value_du_dxx = (unm1(i-1) - 2*unm1(i) + unm1(i+1)) / dx2;
        u(i) = 2.d0*unm1(i) - unm2(i) + C_ijkl(i) * value_du_dxx * (DELTAT^2.d0)/rho(i)+g(i);
    end
    
    
    g(:)=0.d0;
    for i=2:NX
        %UI=ZERO;
    % Corrector 
    %   i) Calcul de f=(A0-K0)*u0 (definition de la source secondaire)
        UI(1:3,1:3,i-1) = [u(i-1) u(i) u(i+1) ;unm1(i-1) unm1(i) unm1(i+1);unm2(i-1) unm2(i) unm2(i+1)];
        %UI=[u(i-1) u(i) u(i+1) ;unm1(i-1) unm1(i) unm1(i+1);unm2(i-1) unm2(i) unm2(i+1)];
    
        g(i)=-sum(sum((dA(1:3,1:3,i-1)-dK(1:3,1:3,i-1)).*UI(1:3,1:3,i-1)))*DELTAT^2/rho(i);
        %g(i)=-sum(sum((dA(1:3,1:3,i-1)-dK(1:3,1:3,i-1)).*UI))*DELTAT^2/rho(ISOURCE);

        
%         unm2(i)=unm1(i);
%         unm1(i)=u(i);
        
    end
    
    %g(:)=0.d0;
    
    
    % Corrector ii) utiliser le nouveau g en tant que la force secondaire
    %               pour calculer u1
    
    for i=2:NX
        value_du_dxx = (u1nm1(i-1) - 2*u1nm1(i) + u1nm1(i+1)) / dx2;
        u1(i) = 2.d0*u1nm1(i) - u1nm2(i) + C_ijkl(i) * value_du_dxx * (DELTAT^2.d0)/rho(i)+g(i);
        %displacement field correction
    end
    
    
    utotal=u+u1; % sommer des u0 et u1 pour ce time step
    % shift des u0
    unm2=unm1;
    unm1=u;
    
    % shift des u1
    u1nm2 = u1nm1;
    u1nm1 = u1;
      
    u=u;
    
    

    
    
    %Dirichlet BC
    u(1) = ZERO;
    u(NX+1) = ZERO;
    
    if SAVE_SEISMOGRAMS
        for irec = 1:NREC
            seisux(it,irec) = u(ix_rec(irec));
        end 
    end
    
    if(mod(it,IT_DISPLAY) == 0)
          if max(u)>umaxv
              umaxv=max(u);
          end   
          plot(xvec,u,'linewidth',2);
          if FE_BOUNDARY
            line([XFE XFE],[-umaxv umaxv],'Color', 'r');
          end

          grid on;
          axis([min(xvec) max(xvec) -umaxv umaxv]);
          xlabel('x, m','fontSize',14);
          ylabel('u','fontSize',14);              
          titlestring = ['TIME STEP = ',num2str(it), ' TIME = ',num2str(t), ' s'];
          title(titlestring ,'fontsize',14);      
          drawnow;

          if PAUSE_ON
              pause(pause_time);
          end
          
          if MAKE_MOVIE
              F_y=getframe(gcf);  %-  capture figure or use gcf to get current figure. Or get current
              writeVideo(vidObj_vy,F_y);  %- add frame to the movie
              %fprintf('Frame for %s captured\n',movie_name_vy);
          end
          
        if DATA_TO_BINARY_FILE && nnz(snapshot_time==it)>0
            filename=[tag 'dis_t_' num2str(it) '.txt'];
            dlmwrite(filename, u);
            fprintf('Data file %s saved to %s\n',filename, pwd);
        end
          %input('Next?');
          
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
  
  if SAVE_SEISMOGRAMS
      for i=1:NREC
          filename=[seis_tag 'x' num2str(xrec(i),'%.2f') 'seis_' num2str(i) '.txt'];
          dlmwrite(filename, [time_vec, seisux(:,i)]);
          fprintf('Seismogram for rec at %.2f saved as %s to %s\n', xrec(i), filename, pwd);
      end
  end
  
  if MAKE_MOVIE
	  close(vidObj_vy);     %- close video file
      fprintf('Video %s saved in %s\n',movie_name_vy, pwd);
  end

  disp('End.');


  
  
  
  
  
  
  
  
  
  
  






