clear all;
close all;

nMrk=6;

dSens=load('input_sensitivity.dat+');
kIndex=dSens(1,1);
phIndex=dSens(2,1);
dkSens=load('input_k.dat+');
ki=dkSens(1,:)';
kiName = strcat('k',num2str(kIndex));

errWt=[1;1;2;2];

lStyles=[': ';':*';'-.';'--';'- '];
colors=['g';'r';'m';'k'];
markers=['^';'o';'s';'v'];

lWidth = [4;4;2;2;2];
mrkSize = 18;

nk = numel(ki);

data=load('../data/SCW-703K.dat');
tData=data(1:5,1);
Data(:,1)=data(1:5,5)/100;
Data(:,2)=data(1:5,3)/100;
Data(:,3)=data(1:5,2)/100;
Data(:,4)=data(1:5,4)/100;
nTimes=numel(tData);

Mtot1=zeros(nTimes,4,nk);

err=zeros(nk,5);

errWtSum = ones(1,4)*errWt;
errWt = errWt/errWtSum;

for ik=1:nk
    ikName = num2str(ik);
    caseName = strcat('case_',kiName,'_',ikName);
    out_dir_name = strcat('results/',caseName);

    fileName = strcat(out_dir_name,'/Time.dat');
    dTime=load(fileName);
    fileName = strcat(out_dir_name,'/Mr.dat');
    Mr=load(fileName);
    fileName = strcat(out_dir_name,'/Mt.dat');
    Mt=load(fileName);
    fileName = strcat(out_dir_name,'/Mck.dat');
    dck=load(fileName);

    t(:,ik) = dTime(:,1);
    Mtot_tmp = Mr + Mt;
    Mtot_tmp(:,end) = dck(:,3);    

    Mtot(:,1,ik)=Mtot_tmp(:,1);
    Mtot(:,2,ik)=Mtot_tmp(:,2) + Mtot_tmp(:,3) + Mtot_tmp(:,4) + Mtot_tmp(:,5) + Mtot_tmp(:,6);
    Mtot(:,3,ik)=Mtot_tmp(:,7) + Mtot_tmp(:,8);
    Mtot(:,4,ik)=Mtot_tmp(:,end);
    
    for i=1:nTimes
        for j=1:numel(t(:,ik))
            if(abs(t(j,ik)-tData(i,1))<1e-3)
                for k=1:4
                    Mtot1(i,k,ik) = Mtot(j,k,ik);
                    err(ik,k) = err(ik,k) + (Mtot1(i,k,ik) - Data(i,k))^2;
                end
                break;
            end
        end
    end

    err(ik,5) = 0;
    for k=1:4
        err(ik,5) = err(ik,5) + errWt(k,1)*err(ik,k);
    end

    fig1 = figure(1);
    hold on;
    for k=1:4
        lStyleColor_tmp = strcat(lStyles(3,:),colors(k,1));
        markerColor_tmp = strcat(markers(k,1),colors(k,1));
        color_tmp = colors(k,1);

        plot(t(:,ik),Mtot(:,k,ik),lStyleColor_tmp,'linewidth',4);
        plot(tData,Mtot1(:,k,ik),markerColor_tmp,'markersize',mrkSize,'markeredgecolor',color_tmp);
        plot(tData,Data(:,k),markerColor_tmp,'markersize',mrkSize, ...
             'markerfacecolor',color_tmp,'markeredgecolor',color_tmp);
    end
    ha=xlabel('time (min)');
    hb=ylabel('mass fraction');
    set([ha hb],'fontsize',30,'fontWeight','bold');
    set(gca,'fontsize',30,'fontWeight','bold');
    NumTicks = 6;
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks));
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks));
    box on;
    filename=sprintf('%s/Mtotal.fig',out_dir_name);
    savefig(fig1,filename);
end

out_dir_name = 'results';
errData = [ki err];
filename=sprintf('%s/err_%s.dat',out_dir_name,kiName);
save(filename,'errData','-ascii');

fig2 = figure(2);
hold on;
for ik=1:nk
    for k=1:4
        lStyleColor_tmp = strcat(lStyles(ik,:),colors(k,1));
        markerColor_tmp = strcat(markers(k,1),colors(k,1));
        color_tmp = colors(k,1);

        plot(t(:,ik),Mtot(:,k,ik),lStyleColor_tmp,'linewidth',lWidth(ik,1));        
        plot(tData,Data(:,k),markerColor_tmp,'markersize',mrkSize, ...
             'markerfacecolor',color_tmp,'markeredgecolor',color_tmp);
    end
end
ha=xlabel('time (min)');
hb=ylabel('mass fraction');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
NumTicks = 6;
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks));
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks));
box on;
filename=sprintf('%s/Mtotal_%s.fig',out_dir_name,kiName);
savefig(fig2,filename);