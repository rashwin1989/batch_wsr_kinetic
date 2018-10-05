clear all;
close all;

%dbstop if naninf

%Species
% 1 - g   -  gas species (C3H8 used)
% 2 - d   -  distillate product (doesn't include gas species)
% 3 - dnr -  distillate non-reactive product
% 4 - map -  MA+ maltene cores (B.P. ~ 480degC)
% 5 - mal -  MA light maltenes (with B.P. < 625degC)
% 6 - mah -  MA heavy maltenes (with B.P. > 625degC)
% 7 - asp -  AS+ asphaltene cores (B.P. ~ 620degC)
% 8 - as  -  AS asphaltenes (taken to be large MW ~ 2000)
% 9 - water/N2
% 10 - ck  - Coke (not part of oil phase)

%%---------------------------------------------------------------------------------%%
%Read and set input parameters
dinp=load('input_main.dat+');
batch=dinp(1,1);              %batch reactor mode
T=dinp(2,1);                  %Temperature in K
P=dinp(3,1);                  %Pressure in Pa
mtot01=dinp(4,1);             %Initial mass of oil phase in g
mtot02=dinp(5,1);             %Initial mass of N2/water phase in g
V0=dinp(6,1);                 %Intial reactor volume in mL
waterInRatio=dinp(7,1);       %Normalized water inlet flow rate
                              %(NWR)
oilInRatio=dinp(8,1);         %Normalized oil inlet flow rate
                              %(NOR)

k_mt=dinp(9,1);               %Two-phase mass transfer co-efficient

tEnd=dinp(10,1);              %End time in mins
tEnd0=dinp(11,1);             %Time scale for water flow rate
                              %normalization in mins
dt=dinp(12,1);                %Time step size in mins
t_write_interval=dinp(13,1);  %Output write interval in mins
err_tol=dinp(14,1);           %Outer loop iterations error
                              %tolerance
chem_tol=dinp(15,1);          %Chemical reaction source term
                              %iterations error tolerance
nIterMax=dinp(16,1);          %Max iterations for outer loop
calcOilIn=dinp(17,1);         %Calculate oil inlet flow rate in code

if(batch==1)
    waterInRatio=0;
    oilInRatio=0;
    nIterMax=1;
end
mdot_in_w=mtot01*waterInRatio/tEnd0; %inlet mass flow rate of water
                                     %in g/min                                     
mdot_in_o = mtot01*oilInRatio/tEnd0; %inlet mass flow rate of oil
                                     %in g/min                  

dmdt_max=[100 100 100 100 100 100 100 100 100]';
initial_equilibrate=0;
out_dir_name = 'results';
mkdir(out_dir_name);
yieldsFile = fopen('results/yields.dat','w');

%%---------------------------------------------------------------------------------%%


%%---------------------------------------------------------------------------------%%
%Read data from files
% chemical kinetic parameters
k1=load('k1.dat+');
k2=load('k2.dat+');
SL1_1=abs(k1(12)); % S_L1 : AS+ - S_L1*(MA+ + MAl + MAh) 
SL2_1=abs(k1(13)); % S_L2 : AS+ - S_L2*(D + Dnr)
SL3_1=abs(k1(14)); % S_L3 : AS+ - S_L3*H2O
SL1_2=abs(k2(12)); % S_L1 : AS+ - S_L1*(MA+ + MAl + MAh) 
SL2_2=abs(k2(13)); % S_L2 : AS+ - S_L2*(D + Dnr)
SL3_2=abs(k2(14)); % S_L3 : AS+ - S_L3*H2O
cgf1=abs(k1(11)); %coke gas fraction: AS+ - cgf*Gas + (1-cgf)*Coke 
cgf2=abs(k2(11)); %coke gas fraction: AS+ - cgf*Gas + (1-cgf)*Coke 

%Thermodynamic parameters
%sp=[21 2 2 8 12 16 14 18 19];
sp=load('species_main.dat+');

d=load('petro.dat+');
Tc=d(sp,4);
Pc=d(sp,5)*1e5;
w =d(sp,6);
MW=d(sp,7);
N=numel(w);
ASindex=N-1;
tk=-sp';

d1=load('y0.dat+');
y1=d1(1,:)';
y2=d1(2,:)';
yoil=y1;
yw=y2;

%%---------------------------------------------------------------------------------%%

for i=1:N
 invMW(i,1) = 1/MW(i,1);    
end
x1=y1.*invMW/(y1'*invMW);
x2=y2.*invMW/(y2'*invMW);
m1=mtot01*y1;
m2=mtot02*y2;
mck1=0; mck2=0; mck=0;
Mtot01 = mtot01;

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
    [dm e_m e_n single_phase] = calc_mass_transfer(1,[100 100 100 100 ...
                        100 100 100 100 100]',m1,m2,P,T,Pc,Tc,w,tk,MW,dt,ASindex);    
    m2 = m2 + dm;
    m1 = m1 - dm;
end
mtot1 = m1'*ones(N,1);
mtot2 = m2'*ones(N,1);
y1 = m1/mtot1;
y2 = m2/mtot2;
x1=y1.*invMW/(y1'*invMW);
x2=y2.*invMW/(y2'*invMW);
mtot = mtot1 + mtot2;
mtot01=mtot1;
mtot02=mtot2;
mtot0=mtot;

rho1=Matlab_density(P,T,Pc,Tc,w,tk,MW,x1)/1000;   %density in g/mL
rho2=Matlab_density(P,T,Pc,Tc,w,tk,MW,x2)/1000;
V1=mtot1/rho1;    %volume in mL
V2=mtot2/rho2;
V0=V1+V2;

% initial guess for outlet mass flow rate
mdot_out = mdot_in_w + mdot_in_o;

mt=zeros(N,1);
m0t=mt;
m01=m1;
m02=m2;
mck01=mck1;
mck02=mck2;
mck0=mck;
m1_prev=m1;
m2_prev=m2;
mdot_out_prev = mdot_out;
y01=y1;
y02=y2;
x01=x1;
x02=x2;
V01=V1;
V02=V2;
rho01=rho1;
rho02=rho2;

Time=0;
Time1=0;
M1=m1';
M2=m2';
Mck1=mck1;
Mck2=mck2;
Mck=mck;
Mr=(m1+m2)';
Mt=mt';
Mtot1=mtot1;
Mtot2=mtot2;
Mdot_out_all=zeros(1,N);
V_oil=V1/(V1+V2);
V_water=V2/(V1+V2);
V_reactor=V1+V2;
reactor_single_phase=single_phase;
R1=zeros(1,N+1);
R2=zeros(1,N+1);
Mdot_intfc=zeros(1,N);
Mdot_out=mdot_out;
Mdot_in_o=mdot_in_o;
Err_mdot=0;
Err_mass=0;
Err_R1=0;
Err_R2=0;
Err_F1=0;
Err_F2=0;
err_cum_mass = 0;
Err_cum_mass = 0;
NIters_outer = 0;
NSteps_chem = 0;
i_timeStep=1;

for t=dt:dt:tEnd

  if (mod(t,t_write_interval)<1e-6)
      disp(['Time: ' num2str(t)]);
  end

  err_iter=1;
  nIter=0;
  single_phase=0;
  %mdot_in_o = -(m1tot_0 - m1tot_00)/dt;  
  mdot_out = mdot_in_w + mdot_in_o;
  mdot_out_prev = mdot_out;
  single_phase_prevIter = single_phase;
  V1_prevIter = V1;
  V2_prevIter = V2;  

  while(err_iter > err_tol & nIter < nIterMax)

      %%---------------------------------------------------------------------------------%%
      
      % advective term contributions
      % if (single_phase)
      %     m2 = m02 + m01 + mdot_in_w*dt*yw + mdot_in_o*dt*yoil - mdot_out*dt*y02;
      %     m1 = m01 - m01;
      %     mt = m0t + mdot_out*dt*y02;
      %     mck1 = mck01;
      %     mck2 = mck02;
      % else
      m2 = m02 + mdot_in_w*dt*yw - mdot_out*dt*y02;
      m1 = m01 + mdot_in_o*dt*yoil;
      mt = m0t + mdot_out*dt*y02;
      mck1 = mck01;
      mck2 = mck02;
          %end

      % prevent -ve species masses
      tmp_m01 = m1;
      tmp_m02 = m2;
      for i=1:N
          if (m1(i,1) < 0)
              m1(i,1) = 0;
          end
          if (m2(i,1) < 0)
              m2(i,1) = 0;
          end
      end
      tmp_m1 = m1;
      tmp_m2 = m2;      
      err_F1 = (tmp_m1 - tmp_m01)';                       
      err_F2 = (tmp_m2 - tmp_m02)';      

      m1tot=m1'*ones(N,1);
      if(m1tot<1e-6)
          single_phase=1;          
      end

      %end of advective term contributions
      %%---------------------------------------------------------------------------------%%

      % chemcial reaction source terms calculation      
      [dm1 dm2 rr1 rr2 err_R1 err_R2 nChemSteps] = calc_chem_source(m1,m2,V1,V2,k1,k2,dt,single_phase,chem_tol);
      m1 = m1 + dm1;
      m2 = m2 + dm2;

      % calculate coke formation in current time step
      masp_ex1 = m1(7,1) - (SL1_1*(m1(4,1)+m1(5,1)+m1(6,1)) + SL2_1*(m1(2,1)+m1(3,1)) + SL3_1*m1(9,1));
      masp_ex2 = m2(7,1) - (SL1_2*(m2(4,1)+m2(5,1)+m2(6,1)) + SL2_2*(m2(2,1)+m2(3,1)) + SL3_2*m2(9,1));
      %masp_ex2 = m2(7,1) - SL3_2*m2(9,1);
      if(masp_ex1 > 0)
          mck1 = mck1 + (1 - cgf1)*masp_ex1;
          m1(1,1) = m1(1,1) + cgf1*masp_ex1;
          m1(7,1) = m1(7,1) - masp_ex1;
          rr1(N+1,1) = (1 - cgf1)*masp_ex1;
          rr1(1,1) = rr1(1,1) + cgf1*masp_ex1;
      end 
      if(masp_ex2 > 0)
          mck2 = mck2 + (1 - cgf2)*masp_ex2;
          m2(1,1) = m2(1,1) + cgf2*masp_ex2;
          m2(7,1) = m2(7,1) - masp_ex2;
          rr2(N+1,1) = (1 - cgf2)*masp_ex2;
          rr2(1,1) = rr2(1,1) + cgf2*masp_ex2;
      end 
      mck = mck1 + mck2;
      
      % some chemical source terms may be > m(i) and -ve
      % some small mass conservation error introduced      
      for i=1:N
          if (m1(i,1) < 0)
              m1(i,1) = 0;
          end
          if (m2(i,1) < 0)
              m2(i,1) = 0;
          end
      end     
      % end of chemical source terms calculation
      %%---------------------------------------------------------------------------------%%       
      % interphase mass transfer calculation
      if (single_phase==0)          

          % calculate interphase mass transfer of species
          [dm e_m e_n single_phase] = calc_mass_transfer(k_mt,dmdt_max,m1,m2,P, ...
                                                         T,Pc,Tc,w,tk,MW,dt,ASindex);
          
          if (single_phase==0)
              m2 = m2 + dm;
              m1 = m1 - dm;                                          
          else
              m2 = m2 + m1;
              dm = m1;
              m1 = m1 - m1;                                          
          end
      else          
          dm = zeros(N,1);
      end
      % interphase mass transfer calculation done
%%---------------------------------------------------------------------------------%%       
      if (single_phase==1)
          mtot2 = m2'*ones(N,1);
          mtot1 = m1'*ones(N,1);
          y2 = m2/mtot2;  
          y1 = y2;
          x2=y2.*invMW/(y2'*invMW);
          x1 = x2;          

          % recalculate mdot_out to conserve reactor volume
          rho2=Matlab_density(P,T,Pc,Tc,w,tk,MW,x2)/1000;
          rho1=rho2;
          V2 = V0;
          V1 = 0;
          if (batch==0)
              mtot2_new = V2*rho2;
              mdot_out = mdot_in_w + mdot_in_o - (mtot1 + mtot2_new + mck - mtot01 - mtot02 - mck0)/dt; 
              %mdot_out = mdot_out + (mtot2 - mtot2_new)/dt;
              
              err_m1 = norm(m1 - m1_prev);              
              err_m2 = norm(m2 - m2_prev);
              err_mdot = norm(mdot_out - mdot_out_prev);
              err_iter = norm([err_m1 err_m2 err_mdot],Inf);

              m1_prev = m1;
              m2_prev = m2;
              mdot_out_prev = mdot_out;
              mtot2 = mtot2_new;
          else              
              V2 = mtot2/rho2;
              err_mdot = 0;
          end
      else
          if(calcOilIn==1)
              mtot1 = m1'*ones(N,1);
              mdot_in_o = mdot_in_o + (mtot01 - mtot1)/dt;
              m1 = m1 + (mtot01 - mtot1)*yoil;
              mtot1 = m1'*ones(N,1);
          else
              mtot1 = m1'*ones(N,1);
          end

          mtot2 = m2'*ones(N,1);

          y1 = m1/mtot1;
          y2 = m2/mtot2;  
          x1=y1.*invMW/(y1'*invMW);
          x2=y2.*invMW/(y2'*invMW);

          % recalculate mdot_out to conserve reactor volume
          rho1=Matlab_density(P,T,Pc,Tc,w,tk,MW,x1)/1000;   %density in g/mL
          rho2=Matlab_density(P,T,Pc,Tc,w,tk,MW,x2)/1000;
          V1 = mtot1/rho1;          
          if (batch==0)
              V2 = V0 - V1;
              mtot2_new = V2*rho2;
              mdot_out = mdot_in_w + mdot_in_o - (mtot1 + mtot2_new + mck - mtot01 - mtot02 - mck0)/dt;               

              err_m1 = norm(m1 - m1_prev);
              err_m2 = norm(m2 - m2_prev);
              err_mdot = norm(mdot_out - mdot_out_prev);
              err_iter = norm([err_m1 err_m2 err_mdot],Inf);

              m1_prev = m1;
              m2_prev = m2;
              mdot_out_prev = mdot_out;
              mtot2 = mtot2_new;
          else
              V2 = mtot2/rho2;
              err_mdot = 0;
          end
      end
      
      nIter = nIter + 1;
      single_phase_prevIter = single_phase;
  end

  mass1 = m1'*ones(N,1);
  mass01 = m01'*ones(N,1);
  mass2 = m2'*ones(N,1);
  mass02 = m02'*ones(N,1);
  err_mass = abs(mass1 + mass2 + mck - mass01 - mass02 - mck0 + mdot_out*dt - ...
                 mdot_in_w*dt - mdot_in_o*dt); 
  err_cum_mass = err_cum_mass + err_mass;%(mass1 + mass2 + mck - mass01 - mass02 - mck0 + mdot_out*dt - ...
                                         %mdot_in_w*dt - mdot_in_o*dt); 

  if (mod(t,t_write_interval)<1e-6)
      Time=[Time;t];
      M1=[M1;m1'];
      M2=[M2;m2'];
      Mr=[Mr;(m1+m2)'];
      Mt=[Mt;mt'];
      Mck1=[Mck1;mck1];
      Mck2=[Mck2;mck2];
      Mck=[Mck;mck];
      Mtot1=[Mtot1;mtot1];
      Mtot2=[Mtot2;mtot2];
      V_oil=[V_oil;V1/(V1+V2)];
      V_water=[V_water;V2/(V1+V2)];
      V_reactor=[V_reactor;V1+V2];
      reactor_single_phase=[reactor_single_phase;single_phase];
      R1=[R1;(rr1/dt)'];
      R2=[R2;(rr2/dt)'];
      Mdot_intfc=[Mdot_intfc;(dm/dt)'];
      Mdot_out=[Mdot_out;mdot_out];
      Mdot_out_all=[Mdot_out_all;mdot_out*y2'];
      Mdot_in_o=[Mdot_in_o;mdot_in_o];
      Err_mass=[Err_mass;err_mass];
      Err_cum_mass=[Err_cum_mass;err_cum_mass];
      Err_mdot=[Err_mdot;err_mdot];          
      Err_R1 = [Err_R1;norm(err_R1,Inf)];
      Err_R2 = [Err_R2;norm(err_R2,Inf)];
      Err_F1 = [Err_F1;norm(err_F1,Inf)];
      Err_F2 = [Err_F2;norm(err_F2,Inf)];
      NIters_outer = [NIters_outer;nIter];
      NSteps_chem = [NSteps_chem;nChemSteps];

      fprintf(['%-10s %-12s %-10s %-6s %-10s %-12s %-6s %-6s %-6s %-6s %-8s %-10s %-10s\n'],'err_mass','err_cum_mass','err_mdot','nIter','nChemSteps','single_phase','mtot1','mtot2','V1','V2','(V1+V2)','mdot_in_o','mdot_out');
      fprintf(['%-10.4e %-12.4e %-10.4e %-6d %-10d %-12d %-6.2f %-6.2f %-6.2f %-6.2f %-8.4e %-10.4e %-10.4e\n'],err_mass,err_cum_mass,err_mdot,nIter,nChemSteps,single_phase,mtot1,mtot2,V1,V2,(V1+V2),mdot_in_o,mdot_out);
  end  

  Mtot01 = Mtot01 + mdot_in_o*dt;
  m01 = m1;
  m02 = m2;
  y02 = y2;
  mck01 = mck1;
  mck02 = mck2;
  mck0 = mck;
  m0t = mt;
  y01 = y1;
  y02 = y2;
  mtot01=mtot1;
  mtot02=mtot2;

  i_timeStep = i_timeStep + 1;
end

%if (mdot_in_o == 0)
M1 = M1/Mtot01;
M2 = M2/Mtot01;
Mck1 = Mck1/Mtot01;
Mck2 = Mck2/Mtot01;
Mck = Mck/Mtot01;
Mt = Mt/Mtot01;
Mr = Mr/Mtot01;
%end

MtFinal = Mt(end,:);
MrFinal = Mr(end,:);
% In MTotFinal, H2O is replaced by Coke as 
% the Nth species
MTotFinal = MtFinal + MrFinal;
MTotFinal(1,end) = Mck(end,1);
% if(mdot_in_o > 0)
%     MTotFinal = MTotFinal/Mtot01;
% end
MTotFinalSum = MTotFinal*ones(N,1);
fprintf('%6.4f\n',MTotFinal);
fprintf('%6.4f\n',MTotFinalSum);
%MTotFinal
%MTotFinalSum

Time_data=[Time Mtot1 Mtot2 Mck1 Mck2 V_oil V_water V_reactor reactor_single_phase Mdot_out Mdot_in_o ...
          Err_mass Err_cum_mass Err_mdot Err_R1 Err_R2 Err_F1 Err_F2 NIters_outer ...
          NSteps_chem];
Mck_data=[Mck1 Mck2 Mck];

fprintf(yieldsFile,'%6.4f\n',MTotFinal);
fprintf(yieldsFile,'%6.4f\n',MTotFinalSum);

filename=sprintf('%s/Time.dat',out_dir_name);
save(filename,'Time_data','-ascii');

filename=sprintf('%s/M1.dat',out_dir_name);
save(filename,'M1','-ascii');
filename=sprintf('%s/M2.dat',out_dir_name);
save(filename,'M2','-ascii');
filename=sprintf('%s/Mck.dat',out_dir_name);
save(filename,'Mck_data','-ascii');

filename=sprintf('%s/Mr.dat',out_dir_name);
save(filename,'Mr','-ascii');
filename=sprintf('%s/Mt.dat',out_dir_name);
save(filename,'Mt','-ascii');

filename=sprintf('%s/Mdot_out_all.dat',out_dir_name);
save(filename,'Mdot_out_all','-ascii');

filename=sprintf('%s/R1.dat',out_dir_name);
save(filename,'R1','-ascii');
filename=sprintf('%s/R2.dat',out_dir_name);
save(filename,'R2','-ascii');

filename=sprintf('%s/Mdot_intfc.dat',out_dir_name);
save(filename,'Mdot_intfc','-ascii');


fclose(yieldsFile);