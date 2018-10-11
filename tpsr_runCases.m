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
T = [683;703];
tEnd = [30;60];
k_mt = [1;0.1;0.05;0.01;0.005;0.001];
NWR = [5;10;25;50;75;100;150;200;300;400];
%Read and set input parameters
dinp=load('input_main.dat+');
batch=dinp(1,1);              %batch reactor mode
%T=dinp(2,1);                  %Temperature in K
P=dinp(3,1);                  %Pressure in Pa
mtot01=dinp(4,1);             %Initial mass of oil phase in g
mtot02=dinp(5,1);             %Initial mass of N2/water phase in g
V0=dinp(6,1);                 %Intial reactor volume in mL
%waterInRatio=dinp(7,1);       %Normalized water inlet flow rate
                              %(NWR)
oilInRatio=dinp(8,1);         %Normalized oil inlet flow rate
                              %(NOR)

%k_mt=dinp(9,1);               %Two-phase mass transfer co-efficient

%tEnd=dinp(10,1);              %End time in mins
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

% if(batch==1)
%     waterInRatio=0;
%     oilInRatio=0;
%     nIterMax=1;
% end
% mdot_in_w=mtot01*waterInRatio/tEnd0; %inlet mass flow rate of water
%                                      %in g/min                                     
% mdot_in_o = mtot01*oilInRatio/tEnd0; %inlet mass flow rate of oil
%                                      %in g/min                  

% dmdt_max=[100 100 100 100 100 100 100 100 100]';
% initial_equilibrate=0;
% out_dir_name = 'results';
% mkdir(out_dir_name);
% yieldsFile = fopen('results/yields.dat','w');

%%---------------------------------------------------------------------------------%%

%Thermodynamic parameters
%sp=[21 2 2 8 12 16 14 18 19];
sp=load('species_main.dat+');

d=load('petro.dat+');
Tc=d(sp,4);
Pc=d(sp,5)*1e5;
w =d(sp,6);
MW=d(sp,7);
% N=numel(w);
% ASindex=N-1;
% tk=-sp';

d1=load('y0.dat+');
y1=d1(1,:)';
y2=d1(2,:)';

%%---------------------------------------------------------------------------------%%
for iT=1:nT

    %Read chemical kinetic parameters data from files
    T_tmp = T(iT,1);
    tEnd_tmp = tEnd(iT,1);
    TName = num2str(T_tmp);
    kFileName = strcat('k1',TName,'.dat+');
    k1=load(kFileName);
    kFileName = strcat('k2',TName,'.dat+');
    k2=load(kFileName);

    for ikmt=1:nkmt

        k_mt_tmp = k_mt(ikmt,1);

        for iNWR=1:nNWR

            NWR_tmp = NWR(iNWR,1);
            iTName = num2str(iT);
            ikmtName = num2str(ikmt);
            iNWRName = num2str(iNWR);
            caseName = strcat('case','_',iTName,'_',ikmtName, ...
                                  '_',iNWRName);
            out_dir_name = strcat('results/',caseName);
            mkdir(out_dir_name);
            yieldsFileName = strcat(out_dir_name,'/yields.dat');
            yieldsFile = fopen(yieldsFileName,'w');

            fprintf(['%-12s %-3s %-4.0f %-5s %-6.3f %-5s %4.0f\n'],caseName,'T=',T_tmp,'kmt=',k_mt_tmp,'NWR=',NWR_tmp);
            
            [MTotFinal MTotFinalSum Time_data M1 M2 Mck_data Mr Mt Mdot_out_all ...
             R1 R2 Mdot_intfc] ... 
            = ...
            twoPhaseStirredReactor ...
            (batch,T_tmp,P,mtot01,mtot02,V0,NWR_tmp,oilInRatio,k_mt_tmp, ...
             tEnd_tmp,tEnd0,dt,t_write_interval,err_tol,chem_tol,nIterMax, ...
             calcOilIn,out_dir_name,k1,k2,sp,Tc,Pc,w,Mw,y1,y2);

            fprintf(['\n']);

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

        end
    end
end
%%---------------------------------------------------------------------------------%%