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
factor = 1.d4;

%Define i of source
xsource=0.01d0*XMAX;    %arbitrary source position in [m]
mindist=XMAX;           
for i=2:NX
    dist=xsource-DELTAX*i;
    if abs(dist)<mindist
        mindist=dist;
        ISOURCE=i;
    end
end
clearvars mindist dist;

fprintf('TIME:\n IT_DISP=%d NSTEP=%d\n', IT_DISPLAY, NSTEP);
fprintf('SPACE:\n NX=%d XMAX=%.2f XMIN=%.2f dX=%.2f\n', NX, XMAX, XMIN, DELTAX);
fprintf(' FE Boundary at XFE=%.2f\n', XFE);
% zero
 ZERO = 0.d0;
 eps=0.00001d0;

% large value for maximum
 HUGEVAL = 1.d+30;

% velocity threshold above which we consider that the code became unstable
 STABILITY_THRESHOLD = 1.d+25;
 
u=zeros(NX+1,1);
u1=zeros(NX+1,1);
xvec=[0:NX]*DELTAX;

unm1=zeros(NX+1,1);
unm2=zeros(NX+1,1);
u1nm1=zeros(NX+1,1);
u1nm2=zeros(NX+1,1);
rho=zeros(NX+1,1);
mu=zeros(NX+1,1);
markers=zeros(NX+1,1);

%find points near the interface
x_trialm1=ZERO;
for i=1:NX+1
    x_trial=(i-1)*DELTAX;
    if (XFE<x_trial) && (XFE>x_trialm1)
        markers(i)=1;
        markers(i-1)=1;
        fprintf('\n FE_BOUNDARY:\n il=%d xl=%.2f < XFE=%.2f < ir=%d xr=%.2f\n', i-1, (i-1)*DELTAX, XFE, i, i*DELTAX);
        eta0=(XFE-(i-2)*DELTAX)/DELTAX;
        eta1=1-eta0;
        fprintf(' eta0=%.2f\n eta1=%.2f\n', eta0, eta1);
    end
    x_trialm1=x_trial;
end

%Define i of source
% xsource=0.01d0*XMAX;    %arbitrary source position in [m]
% mindist=XMAX;           
% for i=2:NX
%     dist=xsource-DELTAX*i;
%     if abs(dist)<mindist
%         mindist=dist;
%         ISOURCE=i;
%     end
% end
% clearvars mindist dist;

% Recievers:
NREC=21;
xdeb=0.2d0*XMAX;
xfin=0.6d0*XMAX;
% %Recievers
% NREC=481;
% xdeb=(i_left-240)*DELTAX;
% xfin=(i_left+240)*DELTAX;
% way=13.50:0.50:53.50;
xrec=zeros(NREC,1);
ix_rec=zeros(NREC,1);
xspacerec = (xfin-xdeb) / double(NREC-1);
for irec=1:NREC
    xrec(irec) = xdeb + double(irec-1)*xspacerec;
end
% [sharedVals,idxsIntoA] = intersect(xrec,way);
% xrec=xrec(idxsIntoA);
% NREC=length(xrec);
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
   fprintf('receiver %d x= %.2f\n',irec,xrec(irec))
   fprintf('closest grid point found at distance %.2f in i = %d\n\n',dist,ix_rec(irec));
end

 
%Apply flat Flat elastic boundary
if FE_BOUNDARY
 fprintf('ON. Lithological boundary\n');
%  cp_above_eb=1800.d0;
%  cp_below_eb=3300.d0;
%  rho_above_eb=2400.d0;
%  rho_below_eb=3200.d0;
 cp_below_eb=1800.d0;
 cp_above_eb=3300.d0;
%  cp_above_eb=5500.d0;
 rho_below_eb=2400.d0;
%  rho_above_eb=3200.d0;
rho_above_eb=2400.d0;
else
 fprintf('OFF. Lithological boundary\n');
%  cp_above_eb=3300.d0;
%  cp_below_eb=3300.d0;
%  rho_above_eb=1800.d0;
%  rho_below_eb=1800.d0;
 cp_above_eb=1800.d0;
 cp_below_eb=1800.d0;
 rho_above_eb=2400.d0;
 rho_below_eb=2400.d0;
end

  %Reflection and transition coefficients
  Refl_coef=(rho_below_eb*cp_below_eb-rho_above_eb*cp_above_eb)/(rho_below_eb*cp_below_eb+rho_above_eb*cp_above_eb);
  Trans_coef=2.d0*rho_below_eb*cp_below_eb/(rho_below_eb*cp_below_eb+rho_above_eb*cp_above_eb);
  if Refl_coef<ZERO
      tmps=', inverse polarity';
  else
      tmps='';
  end
  fprintf('From left to right over the elastic boundary at %.2f m: \n', XFE);
  fprintf('  R= %.2f - reflection%s\n  T= %.2f - transmition\n', Refl_coef, tmps ,Trans_coef);
  clearvars  Refl_coef Trans_coef tmps;

%------------------------------------------------------------------------

  %Initiate video object for vy
  if MAKE_MOVIE
	  movie_name_vy= ['1D_u_' num2str(NX) '_' num2str(DELTAX) '_' num2str(f0) '.avi'];
	  vidObj_vy=VideoWriter(movie_name_vy);
	  open(vidObj_vy);
  end

%------------------------------------------------------------------------

% compute the Lame parameters and density  
% Create Cijkl matrix
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
fprintf('\nStart Cijkl 6D %d elements\n',NX*order^4);
%------------------------------------------------------------------------
for ii = 1:NX+1
    x_trial = (ii-1)*DELTAX;
    if x_trial>XFE
        rho(ii) = rho_above_eb;
        for i=1:order
            for j=1:order
                for k=1:order
                    for l=1:order
                         %C(ii,i,j,k,l)=lambdaa*dkr(i,j)*dkr(k,l)+mua*(dkr(i,k)*dkr(j,l)+dkr(i,l)*dkr(j,k));
                         C(ii,1,1,1,1)=mua;
                         mu(ii)=mua;
                    end
                end
            end
        end
    else
        rho(ii) = rho_below_eb;
        for i=1:order
            for j=1:order
                for k=1:order
                    for l=1:order
                         %C(ii,i,j,k,l)=lambdab*dkr(i,j)*dkr(k,l)+mub*(dkr(i,k)*dkr(j,l)+dkr(i,l)*dkr(j,k));
                         C(ii,1,1,1,1)=mub;
                         mu(ii)=mub;
                    end
                end
            end
        end
    end
end
fprintf('Cijkl(x,y) of size: %s  created\n\n',num2str(size(C)));
 
clearvars densitya cpa csa lambdaa mua densityb cpb csb lambdab mub;


  Courant_number1 = cp_above_eb * DELTAT * DELTAX;
  fprintf('Courant number + = %.4f\n',Courant_number1);
  
  Courant_number2 = cp_below_eb * DELTAT * DELTAX;
  fprintf('Courant number - = %.4f\n\n',Courant_number2);

  if(Courant_number1 > 1.d0) || (Courant_number2 > 1.d0) 
      disp('time step is too large, simulation will be unstable');
      break;
  end

%--------------------------------------------------------------------------
  a = pi*pi*f0*f0;
  umaxv = factor*exp(-a*t0^2)/2.d0*DELTAT;
  dx2 = DELTAX^2.d0;
  v=zeros(size(u));
  total_energy_kinetic = zeros(NSTEP,1);
  total_energy_potential = zeros(NSTEP,1);
  total_energy=zeros(NSTEP,1);
  
  input('Start time loop?')
  for it = 1:NSTEP  
    u(:)=ZERO;
    u1(:)=ZERO;
    g=zeros(NX+1,3);
    g1=zeros(NX+1,3);
    t = double(it-1)*DELTAT;
    % Ricker source time function (second derivative of a Gaussian)
    % source_term = factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);
     source_term = factor * exp(-a*(t-t0)^2);
    % define location of the source
    %i = ISOURCE;
    %u(i) = u(i) + source_term;
    
    
    % Predictor
    
    %unm1(i) = unm1(i) + source_term*DELTAT/rho(i);
    g(ISOURCE)=source_term*DELTAT/rho(ISOURCE);
    
    for i=2:NX
        value_du_dxx = (unm1(i-1) - 2*unm1(i) + unm1(i+1)) / dx2;
        u(i) = 2.d0*unm1(i) - unm2(i) + C(i,1,1,1,1) * value_du_dxx * (DELTAT^2.d0)/rho(i)+g(i);
        
        %v(i) = (u(i) - unm2(i))/(2*DELTAT);
        %value_du_dx = (unm1(i+1) - unm1(i-1)) / (2.d0*DELTAX);
        
        %total_energy_potential(it) = total_energy_potential(it)+0.5*C(i,1,1,1,1) * value_du_dx^2;
    end
   
    g=zeros(NX+1,3);
    g1=zeros(NX+1,3);
    
    A_con = (rho(i)/DELTAT^2)*[0; 1; 0; 0; -2; 0; 0; 1 ; 0];
    K_con = (mu(i)/DELTAX^2)*[0; 0; 0; 1; -2; 1; 0; 0; 0];
    A_opt = (rho(i)/DELTAT^2)*[1/12; 10/12; 1/12; -2/12; -20/12; -2/12; 1/12; 10/12; 1/12];
    K_opt = (mu(i)/DELTAX^2)*[1/12; -2/12; 1/12; 10/12; -20/12; 10/12; 1/12; -2/12; 1/12];
    delta_A = A_opt - A_con;
    delta_K = K_opt - K_con;
    
    for i=2:NX
        %g(i)=f ( (u_0), unm1, unm2 )
        %vs_term =-(A_con-K_con).*stencil;
        %vs_term = sum(vs_term(:,1));
        
        stencil = [0; u(i); 0; unm1(i-1); unm1(i); unm1(i+1); unm2(i-1); unm2(i); unm2(i+1)];
        g = (A_con - K_con).*stencil;
        
    end

    % Corrector - i) calculus of (-dA u_0)
    
     i=ISOURCE ;
     stencilper = [0; u(i); 0; unm1(i-1); unm1(i); unm1(i+1); unm2(i-1); unm2(i); unm2(i+1)];
     g = -(A_opt - K_opt).*stencilper;
    
%    Corrector - ii) calculus of g-modified    

    for i=2:NX
        value_du_dxx = (u1nm1(i-1) - 2*u1nm1(i) + u1nm1(i+1)) / dx2;
        u1(i) = 2.d0*u1nm1(i) - u1nm2(i) + C(i,1,1,1,1) * value_du_dxx * (DELTAT^2.d0)/rho(i)+g1(i);
        
        %v(i) = (u(i) - unm2(i))/(2*DELTAT);
        %value_du_dx = (unm1(i+1) - unm1(i-1)) / (2.d0*DELTAX);
        
        %total_energy_potential(it) = total_energy_potential(it)+0.5*C(i,1,1,1,1) * value_du_dx^2;
    end
    
    for i=2:NX
        stencil = [0; u1(i); 0; u1nm1(i-1); u1nm1(i); u1nm1(i+1); u1nm2(i-1); u1nm2(i); u1nm2(i+1)];
        g1 = -(A_opt-K_opt).*stencil ; 
    end
        
    
    u = u + u1;
    
    
    
    
%     fprintf('it=%d \n umax=%e\n imax=%d\n\n',it,max(u),find(u==max(u)));
    % add the source (force vector located at a given grid point)
    
    total_energy_kinetic(it) = 0.5*sum(rho.*(v.^2));
    
    
    t = double(it-1)*DELTAT;
    % Ricker source time function (second derivative of a Gaussian)
    % source_term = factor * (1.d0 - 2.d0*a*(t-t0)^2)*exp(-a*(t-t0)^2);
    % source_term = factor * exp(-a*(t-t0)^2);
    % define location of the source
    %i = ISOURCE;
    %u(i) = u(i) + source_term;
    %u(i) = u(i) + vs_term;
    
    
    
    
    %Dirichlet BC
    u(1) = ZERO;
    u(NX+1) = ZERO;
    
    if SAVE_SEISMOGRAMS
        for irec = 1:NREC
            seisux(it,irec) = u(ix_rec(irec));
        end 
     end
    
    %Set previous timesteps
    unm2=unm1;
    unm1=u;
    u1nm2=u1nm1;
    u1nm1=u1;
    
    
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
