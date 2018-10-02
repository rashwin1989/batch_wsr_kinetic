clear all;
close all;

%%---------------------------------------------------------------------------------%%
%Set input parameters
batch=0;
T=703;    %temperature in K
P=30e06;  %pressure in Pa
mtot01=10; %mass in g
mtot02=25;
%Mtot01=mtot01;
tEnd=15; %time in mins
tEnd0=60;
dt=0.01;
t_write_interval=0.01;
err_tol=1e-09;
chem_tol=1e-12;
nIterMax=108;
i_debug_timeStep=770;
ck0_factor = 1e-12;
ck_phase_model = 3;
V0=100;
waterInRatio=100;
if(batch==1)
    waterInRatio=0;
    nIterMax=1;
end
mdot_in_w=mtot01*waterInRatio/tEnd0; %inlet mass flow rate of water
                                 %g/min
                                 %mdot_in_w = 24/tEnd0;
mdot_in_o = 0/tEnd0;
Mtot01 = mdot_in_o*tEnd + mtot01; 
k_mt=1;
dmdt_max=[0.01 100 100 100]';
initial_equilibrate=1;
out_dir_name = 'results';
mkdir(out_dir_name);
logfileID = fopen('logwsr','w');

%%---------------------------------------------------------------------------------%%


%%---------------------------------------------------------------------------------%%
%Read data from files
% chemical kinetic parameters
k=load('k703.dat'); 

%Thermodynamic parameters
sp=[6 15 18 19];

d=load('petro.dat+.C3H8');
Tc=d(sp,4);
Pc=d(sp,5)*1e5;
w =d(sp,6);
MW=d(sp,7);
N=numel(MW);
tk=-sp';

y1=load('y01.dat');
y1=y1';
yoil=zeros(N+3,1);
for i=1:N
    yoil(i,1) = y1(i,1);
end
y2=[zeros(N-1,1);1];

%%---------------------------------------------------------------------------------%%

for i=1:N
 invMW(i,1) = 1/MW(i,1);    
end
x1=y1.*invMW/(y1'*invMW);
x2=y2.*invMW/(y2'*invMW);
m1=mtot01*y1;
m2=mtot02*y2;
%1:maltene_dist; 2:maltene_vr; 3:asphaltene; 4:h2o

% determine initial volume of water phase in reactor
% oil phase and water phase with input compositions
rho01=Matlab_density(P,T,Pc,Tc,w,tk,MW,x1)/1000;   %density in g/mL
rho02=Matlab_density(P,T,Pc,Tc,w,tk,MW,x2)/1000;
V01=mtot01/rho01;
V02=mtot02/rho02;
V0=V01+V02;
mtot0=mtot01+mtot02;

% determine initial compositions of oil and water phase using phase
% equilibrium calculation
single_phase = 0;
if (initial_equilibrate==1)
    [dm e_m e_n single_phase] = calc_mass_transfer(1,[100 100 100 ...
                        100]',m1,m2,P,T,Pc,Tc,w,tk,MW,dt);    
    m2 = m2 + dm;
    m1 = m1 - dm;
    if (abs(dm(1,1))>dmdt_max(1,1)*dt)
        m2(1,1) = m2(1,1) - dm(1,1);
        m1(1,1) = m1(1,1) + dm(1,1);
        dm(1,1) = dm(1,1)/abs(dm(1,1))*dmdt_max(1,1)*dt;
        m2(1,1) = m2(1,1) + dm(1,1);
        m1(1,1) = m1(1,1) - dm(1,1);
    end
end
mtot1 = m1'*ones(N,1);
mtot2 = m2'*ones(N,1);
y1 = m1/mtot1;
y2 = m2/mtot2;
x1=y1.*invMW/(y1'*invMW);
x2=y2.*invMW/(y2'*invMW);
mtot = mtot1 + mtot2;

rho1=Matlab_density(P,T,Pc,Tc,w,tk,MW,x1)/1000;   %density in g/mL
rho2=Matlab_density(P,T,Pc,Tc,w,tk,MW,x2)/1000;
V1=mtot1/rho1;    %volume in mL
V2=mtot2/rho2;
V0=V1+V2;

% arrays to store masses of all species including minor species 
% Gas, coke, Maltene_dist_nr
mc1=zeros(N+3,1);
mc2=zeros(N+3,1);
mct=zeros(N+3,1);
%5:maltene_dist_nr; 6:gas; 7:coke
for i=1:numel(y1)
    mc1(i,1)=mtot1*y1(i,1);
    mc2(i,1)=mtot2*y2(i,1);
end

%setting initial coke mass in phases
mck0=mtot1*ck0_factor;
gamma=k(11);
%ratio=gamma/V2/(1/V1 + gamma/V2);
if(ck_phase_model==1)
    ratio=V1/(V1+V2)*gamma;
else
    if (ck_phase_model==3)
        ratio=gamma/V2/(1/(V1 + 1e-12) + gamma/V2);
    end
    ratio=mc1(3,1)/(mc1(3,1)+mc2(3,1));
end
if (ratio>1) ratio=1; end
mc1(7,1) = mck0*ratio;
mc2(7,1) = mck0 - mc1(7,1);

mtot1 = mc1'*ones(N+3,1);
mtot2 = mc2'*ones(N+3,1);
yc1 = mc1/mtot1;
yc2 = mc2/mtot2;

% initial guess for outlet mass flow rate
mdot_out = mdot_in_w + mdot_in_o;

m01=m1;
m02=m2;
mc01=mc1;
mc02=mc2;
mc0t=mct;
mc1_prev=mc1;
mc2_prev=mc2;
mdot_out_prev = mdot_out;
mc0t=mct;
y01=y1;
y02=y1;
yc01=yc1;
yc02=yc2;
x01=x1;
x02=x2;
V01=V1;
V02=V2;
rho01=rho1;
rho02=rho2;

yw=zeros(N+3,1);
yw(4,1)=1;

Time=0;
Time1=0;
Mc1=mc1;
Mc2=mc2;
Mr=mc1+mc2;
Mt=mct;
V_oil=V1/(V1+V2);
V_water=V2/(V1+V2);
V_reactor=V1+V2;
reactor_single_phase=single_phase;
R1=zeros(N+3,1);
R2=zeros(N+3,1);
Mdot_intfc=zeros(N+3,1);
Mdot_out=mdot_out;
Err_mdot=0;
Err_mass=0;
Err_R1=zeros(N+3,1);
Err_R2=zeros(N+3,1); 
Err_F1=zeros(N+3,1);
Err_F2=zeros(N+3,1); 
Err_cum_mass = 0;
NIters_outer = 0;
NSteps_chem = 0;
i_timeStep=1;

for t=dt:dt:tEnd

  if (mod(t,t_write_interval)<1e-6)
      display(['Time: ' num2str(t)]);
  end
  err_iter=1;
  nIter=0;  
  single_phase=0;
  mdot_out = mdot_in_w + mdot_in_o;

  while(err_iter > err_tol & nIter < nIterMax)

      if (i_timeStep == i_debug_timeStep)
          fprintf(logfileID, '%5d\r\n\n', nIter);
          fprintf(logfileID, '%12.8f\r\n\n', mc01);
          fprintf(logfileID, '%12.8f\r\n\n', mc02);
          fprintf(logfileID, '%12.8f\r\n\n', mdot_in_w);
          fprintf(logfileID, '%12.8f\r\n\n', mdot_out);
          fprintf(logfileID, '%12.8f\r\n\n', yc2);
      end

      mc2 = mc02 + mdot_in_w*dt*yw - mdot_out*dt*yc2;
      mc1 = mc01 + mdot_in_o*dt*yoil;
      mct = mc0t + mdot_out*dt*yc2;

      if (i_timeStep == i_debug_timeStep)
          fprintf(logfileID, '%12.8f\r\n\n', mc1);
          fprintf(logfileID, '%12.8f\r\n\n', mc2);
      end
      
      % prevent -ve species masses
      tmp_mc01 = mc1;
      tmp_mc02 = mc2;
      for i=1:(N+3)
          if (mc1(i,1) < 0)
              mc1(i,1) = 0;
          end
          if (mc2(i,1) < 0)
              mc2(i,1) = 0;
          end
      end
      tmp_mc1 = mc1;
      tmp_mc2 = mc2;      
      err_F1 = tmp_mc1 - tmp_mc01;                       
      err_F2 = tmp_mc2 - tmp_mc02;      

      mc1tot=mc1'*ones(N+3,1);
      if((mc1tot-mc1(7,1))<1e-6)
          single_phase=1;
          mc2(7,1) = mc2(7,1) + mc1(7,1);
          mc1(7,1) = mc1(7,1) - mc1(7,1);
      end

      dtTotal = 0;
      dt1 = dt;
      nChemSteps = 0;
      rr1 = zeros(N+3,1);
      rr2 = zeros(N+3,1);
      while (dtTotal < dt)
          if (single_phase==1)
              [dm2] = calc_chem_RR(V2,mc2,k,dt1);
              dm1 = zeros(numel(mc1),1);
          else
              if (V1<1e-8)
                  for i=1:N
                      m1(i,1) = mc1(i,1);
                  end 
                  m1(1,1) = m1(1,1) + mc1(5,1) + mc1(6,1);
                  mtot1 = m1'*ones(N,1);
                  y1 = m1/mtot1;
                  x1=y1.*invMW/(y1'*invMW);
                  rho1=Matlab_density(P,T,Pc,Tc,w,tk,MW,x1)/1000;
                  mtot1 = mc1'*ones(N+3,1);
                  V1=mtot1/rho1;
              end
              [dm1] = calc_chem_RR(V1,mc1,k,dt1);
              [dm2] = calc_chem_RR(V2,mc2,k,dt1);
          end
          dtFactorMin = 1;
          for i=1:numel(mc2)
              if ((mc2(i,1) + dm2(i,1)) < -chem_tol)
                  dtFactor = mc2(i,1)/(abs(dm2(i,1)) + 1e-12);
                  dtFactorMin = min(dtFactorMin,dtFactor);
              end
          end
          dt1 = dt1*dtFactorMin;
          dm1 = dm1*dtFactorMin;
          dm2 = dm2*dtFactorMin;
          mc1 = mc1 + dm1;
          mc2 = mc2 + dm2;
          rr1 = rr1 + dm1;
          rr2 = rr2 + dm2;

          tmp_mc01 = mc1;
          tmp_mc02 = mc2;
          for i=1:(N+3)
              if (mc1(i,1) < 0)
                  mc1(i,1) = 0;
              end
              if (mc2(i,1) < 0)
                  mc2(i,1) = 0;
              end
          end
          tmp_mc1 = mc1;
          tmp_mc2 = mc2;
          err_R1 = tmp_mc1 - tmp_mc01;
          err_R2 = tmp_mc2 - tmp_mc02;

          dtTotal = dtTotal + dt1;
          dt1 = dt - dtTotal;
          nChemSteps = nChemSteps + 1;
      end
      
      % some chemical source terms may be > m(i) and -ve
      % some small mass conservation error introduced      
      for i=1:(N+3)
          if (mc1(i,1) < 0)
              mc1(i,1) = 0;
          end
          if (mc2(i,1) < 0)
              mc2(i,1) = 0;
          end
      end      
      
      if (single_phase==0)
          % construct masses of 4 species for phase equilibrium
          % species 1: maltene_dist combines ma_dr, ma_dnr, gas
          for i=1:N
              m1(i,1) = mc1(i,1);
              m2(i,1) = mc2(i,1);
          end
          m1(1,1) = m1(1,1) + mc1(5,1) + mc1(6,1);
          may1(1,1) = mc1(1,1)/(m1(1,1)+1e-12);
          may1(2,1) = mc1(5,1)/(m1(1,1)+1e-12);
          may1(3,1) = mc1(6,1)/(m1(1,1)+1e-12);
          m2(1,1) = m2(1,1) + mc2(5,1) + mc2(6,1);
          may2(1,1) = mc2(1,1)/(m2(1,1)+1e-12);
          may2(2,1) = mc2(5,1)/(m2(1,1)+1e-12);
          may2(3,1) = mc2(6,1)/(m2(1,1)+1e-12);
          
          % coke redistribution according to gamma ratio
          mcktot = mc1(7,1) + mc2(7,1);
          %ratio=gamma/V2/(1/V1 + gamma/V2);
          if(ck_phase_model==1)
              ratio=V1/(V1+V2)*gamma;
          else
              if (ck_phase_model==3)
                  ratio=gamma/V2/(1/(V1+1e-12) + gamma/V2);
              end
              ratio=mc1(3,1)/(mc1(3,1)+mc2(3,1));
          end          
          if (ratio>1) ratio=1; end
          mc1(7,1) = mcktot*ratio;
          mc2(7,1) = mcktot - mc1(7,1);

          % calculate interphase mass transfer of species
          [dm e_m e_n single_phase] = calc_mass_transfer(k_mt,dmdt_max,m1,m2,P, ...
                                                         T,Pc,Tc,w,tk,MW,dt);
          
          if (single_phase==0)
              m2 = m2 + dm;
              m1 = m1 - dm;

              if (abs(dm(1,1))>dmdt_max(1,1)*dt)
                  m2(1,1) = m2(1,1) - dm(1,1);
                  m1(1,1) = m1(1,1) + dm(1,1);
                  dm(1,1) = dm(1,1)/abs(dm(1,1))*dmdt_max(1,1)*dt;
                  m2(1,1) = m2(1,1) + dm(1,1);
                  m1(1,1) = m1(1,1) - dm(1,1);
              end
              
              dmc = zeros(N+3,1);
              for i=1:N
                  dmc(i,1) = dm(i,1);
              end
              if(dm(1,1)>=0)
                  dmc(1,1) = may1(1,1)*dm(1,1);
                  dmc(5,1) = may1(2,1)*dm(1,1);
                  dmc(6,1) = may1(3,1)*dm(1,1);
              else
                  dmc(1,1) = may2(1,1)*dm(1,1);
                  dmc(5,1) = may2(2,1)*dm(1,1);
                  dmc(6,1) = may2(3,1)*dm(1,1);
              end
          else
              mc2 = mc2 + mc1;
              mc1 = mc1 - mc1;
              for i=1:N
                  m1(i,1) = 0;
                  m2(i,1) = mc2(i,1);
              end          
              m2(1,1) = m2(1,1) + mc2(5,1) + mc2(6,1);          
              dmc = zeros(N+3,1);
          end
      else
          for i=1:N
              m1(i,1) = 0;
              m2(i,1) = mc2(i,1);
          end          
          m2(1,1) = m2(1,1) + mc2(5,1) + mc2(6,1);          
          dmc = zeros(N+3,1);
      end

      if (single_phase==1)
          mtot2 = m2'*ones(N,1);
          mtot1 = m1'*ones(N,1);
          y2 = m2/mtot2;  
          y1 = y2;
          x2=y2.*invMW/(y2'*invMW);
          x1 = x2;
          mc2 = mc2 + mc1;
          mc1 = mc1 - mc1;
          mtot2 = mc2'*ones(N+3,1);
          mtot1 = mc1'*ones(N+3,1);
          yc2 = mc2/mtot2;
          yc1 = yc2;

          % recalculate mdot_out to conserve reactor volume
          rho2=Matlab_density(P,T,Pc,Tc,w,tk,MW,x2)/1000;
          rho1=rho2;
          V2 = V0;
          V1 = 0;
          if (batch==0)
              mtot2_new = V2*rho2;
              mdot_out = mdot_out + (mtot2 - mtot2_new)/dt;
              
              err_mc1 = norm(mc1 - mc1_prev);              
              err_mc2 = norm(mc2 - mc2_prev);
              err_mdot = norm(mdot_out - mdot_out_prev);
              err_iter = norm([err_mc1 err_mc2 err_mdot],Inf);

              mc1_prev = mc1;
              mc2_prev = mc2;
              mdot_out_prev = mdot_out;
          else              
              V2 = mtot2/rho2;
          end
      else
          mtot1 = m1'*ones(N,1);
          mtot2 = m2'*ones(N,1);
          y1 = m1/mtot1;
          y2 = m2/mtot2;  
          x1=y1.*invMW/(y1'*invMW);
          x2=y2.*invMW/(y2'*invMW);
 
          mc1 = mc1 - dmc;
          mc2 = mc2 + dmc;

          mtot1 = mc1'*ones(N+3,1);
          mtot2 = mc2'*ones(N+3,1);
          yc1 = mc1/mtot1;
          yc2 = mc2/mtot2;

          % recalculate mdot_out to conserve reactor volume
          rho1=Matlab_density(P,T,Pc,Tc,w,tk,MW,x1)/1000;   %density in g/mL
          rho2=Matlab_density(P,T,Pc,Tc,w,tk,MW,x2)/1000;
          V1 = mtot1/rho1;          
          if (batch==0)
              V2 = V0 - V1;
              mtot2_new = V2*rho2;
              mdot_out = mdot_out + (mtot2 - mtot2_new)/dt;

              err_mc1 = norm(mc1 - mc1_prev);
              err_mc2 = norm(mc2 - mc2_prev);
              err_mdot = norm(mdot_out - mdot_out_prev);
              err_iter = norm([err_mc1 err_mc2 err_mdot],Inf);

              mc1_prev = mc1;
              mc2_prev = mc2;
              mdot_out_prev = mdot_out;
          else
              V2 = mtot2/rho2;
          end
      end

      nIter = nIter + 1;

  end

  mass1 = mc1'*ones(N+3,1);
  mass01 = mc01'*ones(N+3,1);
  mass2 = mc2'*ones(N+3,1);
  mass02 = mc02'*ones(N+3,1);
  err_mass = abs(mass1 + mass2 - mass01 - mass02 + mdot_out*dt - ...
                 mdot_in_w*dt - mdot_in_o*dt); 
  Err_cum_mass = Err_cum_mass + (mass1 + mass2 - mass01 - mass02 + mdot_out*dt - ...
                 mdot_in_w*dt - mdot_in_o*dt);
 
  if (mod(t,t_write_interval)<1e-6)
      Time=[Time t];
      Mc1=[Mc1 mc1];
      Mc2=[Mc2 mc2];
      Mr=[Mr mc1+mc2];
      Mt=[Mt mct];
      V_oil=[V_oil V1/(V1+V2)];
      V_water=[V_water V2/(V1+V2)];
      V_reactor=[V_reactor V1+V2];
      reactor_single_phase=[reactor_single_phase single_phase];
      R1=[R1 rr1/dt];
      R2=[R2 rr2/dt];
      Mdot_intfc=[Mdot_intfc dmc/dt];
      Mdot_out=[Mdot_out mdot_out];           
  end
  Time1=[Time1 t];
  Err_mass=[Err_mass err_mass];
  if (batch==0)
      Err_mdot=[Err_mdot err_mdot];          
  end 
  Err_R1 = [Err_R1 err_R1];
  Err_R2 = [Err_R2 err_R2];
  Err_F1 = [Err_F1 err_F1];
  Err_F2 = [Err_F2 err_F2];
  NIters_outer = [NIters_outer nIter];
  NSteps_chem = [NSteps_chem nChemSteps];

  mc01 = mc1;
  mc02 = mc2;
  mc0t = mct;
  yc01 = yc1;
  yc02 = yc2;
  mtot01=mtot1;
  mtot02=mtot2;

  i_timeStep = i_timeStep + 1;
end

Mc1 = Mc1/Mtot01;
Mc2 = Mc2/Mtot01;
Mt = Mt/Mtot01;
Mr = Mr/Mtot01;
MtFinal = Mt(:,end);
MrFinal = Mr(:,end);
MtFinal(1,1)=MtFinal(1,1) + Mc2(1,end);
MtFinal(5,1)=MtFinal(5,1) + Mc2(5,end);
MtFinal(6,1)=MtFinal(6,1) + Mc2(6,end);
MrFinal(1,1)=MrFinal(1,1) - Mc2(1,end);
MrFinal(5,1)=MrFinal(5,1) - Mc2(5,end);
MrFinal(6,1)=MrFinal(6,1) - Mc2(6,end);
MTotFinal = MtFinal + MrFinal;


filename=sprintf('%s/Time.dat',out_dir_name);
save(filename,'Time','-ascii');
filename=sprintf('%s/Mc1.dat',out_dir_name);
save(filename,'Mc1','-ascii');
filename=sprintf('%s/Mc2.dat',out_dir_name);
save(filename,'Mc2','-ascii');
filename=sprintf('%s/Mr.dat',out_dir_name);
save(filename,'Mr','-ascii');
filename=sprintf('%s/Mt.dat',out_dir_name);
save(filename,'Mt','-ascii');
filename=sprintf('%s/R1.dat',out_dir_name);
save(filename,'R1','-ascii');
filename=sprintf('%s/R2.dat',out_dir_name);
save(filename,'R2','-ascii');
filename=sprintf('%s/Mdot_intfc.dat',out_dir_name);
save(filename,'Mdot_intfc','-ascii');
filename=sprintf('%s/Mdot_out.dat',out_dir_name);
save(filename,'Mdot_out','-ascii');
filename=sprintf('%s/Err_mdot.dat',out_dir_name);
save(filename,'Err_mdot','-ascii');
filename=sprintf('%s/Err_mass.dat',out_dir_name);
save(filename,'Err_mass','-ascii');


fig1 = figure(1);          
hold on;
plot(Time(1,:), Mt(1,:)+Mt(5,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
plot(Time(1,:), Mt(2,:), '-rs' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'r');
plot(Time(1,:), Mt(3,:), '-md' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'm');
plot(Time(1,:), Mt(7,:), '-kv' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'k');
plot(Time(1,:), Mt(6,:), '-c^' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'c');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('mass fraction');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/M_tank.fig',out_dir_name);
savefig(fig1,filename);


fig2 = figure(2);          
hold on;
plot(Time(1,:), Mr(1,:)+Mr(5,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
plot(Time(1,:), Mr(2,:), '-rs' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'r');
plot(Time(1,:), Mr(3,:), '-md' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'm');
plot(Time(1,:), Mr(7,:), '-kv' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'k');
plot(Time(1,:), Mr(6,:), '-c^' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'c');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('mass fraction');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/M_reactor.fig',out_dir_name);
savefig(fig2,filename);


fig3 = figure(3);          
hold on;
plot(Time(1,:), Mc1(1,:)+Mc1(5,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
plot(Time(1,:), Mc1(2,:), '-rs' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'r');
plot(Time(1,:), Mc1(3,:), '-md' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'm');
plot(Time(1,:), Mc1(7,:), '-kv' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'k');
plot(Time(1,:), Mc1(6,:), '-c^' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'c');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('mass fraction');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/M_reactor_oil.fig',out_dir_name);
savefig(fig3,filename);


fig4 = figure(4);          
hold on;
plot(Time(1,:), Mc2(1,:)+Mc2(5,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
plot(Time(1,:), Mc2(2,:), '-rs' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'r');
plot(Time(1,:), Mc2(3,:), '-md' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'm');
plot(Time(1,:), Mc2(7,:), '-kv' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'k');
plot(Time(1,:), Mc2(6,:), '-c^' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'c');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('mass fraction');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/M_reactor_water.fig',out_dir_name);
savefig(fig4,filename);


fig5 = figure(5);
hold on;
plot(Time(1,:), V_oil(1,:), '-rs' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'r');
plot(Time(1,:), V_water(1,:), '-bo' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'b');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('volume fraction');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/V_phases.fig',out_dir_name);
savefig(fig5,filename);

fig6 = figure(6);
hold on;
plot(Time(1,:), V_reactor(1,:), '-ks' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'k');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('volume (mL)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/V_reactor.fig',out_dir_name);
savefig(fig6,filename);


fig7 = figure(7);
hold on;
barY = [MtFinal(1,1)+MtFinal(5,1) MtFinal(2,1) MtFinal(3,1) MtFinal(6,1) ...
        MtFinal(7,1); MrFinal(1,1)+MrFinal(5,1) MrFinal(2,1) ...
        MrFinal(3,1) MrFinal(6,1) MrFinal(7,1)];
h=bar(barY, 'stacked', 'lineWidth', 2);
set(h(1), 'facecolor', 'y'); 
set(h(2), 'facecolor', 'r');
set(h(3), 'facecolor', 'm');
set(h(4), 'facecolor', 'c');
set(h(5), 'facecolor', 'k');
hb=ylabel('mass fraction');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
NumTicks = 2;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks));
set(gca, 'XTickLabel', {'tank','reactor'});
box on;
filename=sprintf('%s/Product_yields.fig',out_dir_name);
savefig(fig7,filename);


fig8 = figure(8);          
hold on;
plot(Time(1,:), R1(1,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
plot(Time(1,:), R1(5,:), '-yo' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'y');
plot(Time(1,:), R1(2,:), '-rs' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'r');
plot(Time(1,:), R1(3,:), '-md' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'm');
plot(Time(1,:), R1(7,:), '-kv' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'k');
plot(Time(1,:), R1(6,:), '-c^' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'c');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('reaction rate (g/min)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/R1.fig',out_dir_name);
savefig(fig8,filename);

fig9 = figure(9);          
hold on;
plot(Time(1,:), R2(1,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
plot(Time(1,:), R2(5,:), '-yo' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'y');
plot(Time(1,:), R2(2,:), '-rs' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'r');
plot(Time(1,:), R2(3,:), '-md' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'm');
plot(Time(1,:), R2(7,:), '-kv' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'k');
plot(Time(1,:), R2(6,:), '-c^' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'c');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('reaction rate (g/min)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/R2.fig',out_dir_name);
savefig(fig9,filename);

R = R1 + R2;
fig10 = figure(10);          
hold on;
plot(Time(1,:), R(1,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
plot(Time(1,:), R(5,:), '-yo' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'y');
plot(Time(1,:), R(2,:), '-rs' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'r');
plot(Time(1,:), R(3,:), '-md' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'm');
plot(Time(1,:), R(7,:), '-kv' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'k');
plot(Time(1,:), R(6,:), '-c^' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'c');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('reaction rate (g/min)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/R.fig',out_dir_name);
savefig(fig10,filename);


fig11 = figure(11);          
hold on;
plot(Time(1,:), Mdot_intfc(1,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
plot(Time(1,:), Mdot_intfc(5,:), '-yo' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'y');
plot(Time(1,:), Mdot_intfc(2,:), '-rs' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'r');
plot(Time(1,:), Mdot_intfc(3,:), '-md' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'm');
plot(Time(1,:), Mdot_intfc(7,:), '-kv' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'k');
plot(Time(1,:), Mdot_intfc(6,:), '-c^' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'c');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('mass flux (g/min)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/Mdot_intfc.fig',out_dir_name);
savefig(fig11,filename);

fig12 = figure(12);          
hold on;
plot(Time(1,:), Mdot_out(1,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('mass flow rate (g/min)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/Mdot_out.fig',out_dir_name);
savefig(fig12,filename); 

fig13 = figure(13);          
hold on;
plot(Time1(1,:), Err_mdot(1,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('mass error (g)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/Err_mdot.fig',out_dir_name);
savefig(fig13,filename); 

fig14 = figure(14);          
hold on;
plot(Time1(1,:), Err_mass(1,:), '-go' ,  'linewidth', 2, ...
     'markersize', 6, 'markerfacecolor', 'g');
xlim([Time(1,1) Time(1,end)]);
ha=xlabel('time (min)');
hb=ylabel('mass error (g)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
box on;
filename=sprintf('%s/Err_mass.fig',out_dir_name);
savefig(fig14,filename); 


fclose(logfileID);