%%%%%_______________PML(Perfectly Matched Layer)_using_FDTD__________%%%%%%%
%_________________________SOUGATA_CHATTERJEE________________________________
%%%%%%%%%_________________SAMEER_KOLKATA_CENTER___________________%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%________02/03/2011_________%%%%%%%%%%%%%
ie=60;
je=60;
ic=fix(ie/2);
jc=fix(je/2);
epsz=8.8e-12;
ddx=0.1;
dt=ddx/6e8;
%%%%%%-----Time_interval-----------%%%%%%%%%%%%
nsteps=100;
%%%%%%-----Time_interval-----------%%%%%%%%%%%%
ez=zeros(ie,je);
dz=zeros(ie,je);
hy=zeros(ie,je);
ihx=zeros(ie,je);
ihy=zeros(ie,je);
hx=zeros(ie,je);
ga=ones(ie,je);
%%%%%%-----wave_specification-----------%%%%%%%%%%%%
t0=40.0;
spread=15.0;
%------------calculated the PML parameter---------------

    gi2=ones(ie);
    gi3=ones(ie);
    fi1=zeros(ie);
    fi2=ones(ie);
    fi3=ones(ie);

    gj2=ones(ie);
    gj3=ones(ie);
    fj1=zeros(ie);
    fj2=ones(ie);
    fj3=ones(ie);
    
npml=20;
for i=1:npml
    xnum=npml-i;
    xd=npml;
    xxn=xnum/xd;
    xn=0.33*((xxn)^3);
    gj2(i)=1.0/(1.0+xn);
    gi2(ie-1-i)=1.0/(1.0+xn);
    gi3(i)=(1.0-xn)/(1.0+xn);
    gi3(ie-i-1)=(1.0-xn)/(1.0+xn);
    xxn=(xnum-0.5)/xd;
    xn=0.25*((xxn)^3);
    fi1(i)=xn;
    fi1(ie-2-i)=xn;
    fi2(i)=1.0/(1.0+xn);
    fi2(ie-2-i)=1.0/(1.0+xn);
    fi3(i)=(1.0-xn)/(1.0+xn);
    fi3(ie-2-i)=(1.0-xn)/(1.0+xn);
end
for j=1:npml
    xnum=npml-j;
    xd=npml;
    xxn=xnum/xd;
    xn=0.33*((xxn)^3);
    gj2(j)=1.0/(1.0+xn);
    gi2(je-1-j)=1.0/(1.0+xn);
    gj3(j)=(1.0-xn)/(1.0+xn);
    gj3(je-j-1)=(1.0-xn)/(1.0+xn);
    xxn=(xnum-0.5)/xd;
    xn=0.25*((xxn)^3);
    fj1(j)=xn;
    fj1(je-2-j)=xn;
    fj2(j)=1.0/(1.0+xn);
    fj2(je-2-j)=1.0/(1.0+xn);
    fj3(j)=(1.0-xn)/(1.0+xn);
    fj3(je-2-j)=(1.0-xn)/(1.0+xn);
end
T=0;
%%%%%%------------Main Loop Begins----------------%%%%%%%%%%%%
for n=1:nsteps %ok<ALIGN>
    T=T+1;
    
    for j=2:je
        for i=2:ie
            dz(i,j)=gi3(i)*gj3(j)*dz(i,j)+gi2(i)*gj2(j)*0.5*(hy(i,j)-hy(i-1,j)-hx(i,j)+hx(i,j-1));
        end
    end
   
     pulse=-2.0*((t0-T)./spread).*exp(-1.*((t0-T)./spread)^2);
%     pulse=exp(-1.*((t0-T)./spread)^2);
    dz(ic,jc)=dz(ic,jc)+pulse;
    for j=2:je
        for i=2:ie
            ez(i,j)=ga(i,j).*dz(i,j);
        end
    end
%%%%%%-----ABC-----------%%%%%%%%%%%%
  for j=1:je-1
     ez(1,j)=0;
     ez(ie-1,j)=0;
  end
 for i=1:ie-1
     ez(i,1)=0;
     ez(i,je-1)=0;
 end     
 %%%%%%-----ABC-----------%%%%%%%%%%%%
    for j=1:je-1
        for i=1:ie-1
            curl_e=ez(i,j)-ez(i,j+1);
            ihx(i,j)=ihx(i,j)+fi1(i)*curl_e;
            hx(i,j)=fj3(j)*hx(i,j)+fj2(j)*0.5*(curl_e+ihx(i,j));
        end
    end
       
     for j=1:je-1
        for i=1:ie-1
            curl_e=ez(i+1,j)-ez(i,j);
            ihy(i,j)=ihy(i,j)+fj1(j)*curl_e;
            hy(i,j)=fi3(i)*hy(i,j)+fi2(i)*0.5*(curl_e+ihy(i,j));
        end
    end
 imagesc(ez);
 title(['Time = ',num2str(n)]);
 colorbar
 pause(0.02);
end