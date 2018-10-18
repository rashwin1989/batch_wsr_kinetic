clear all;
close all;

nMrk=6;

dinp=load('plot_input.dat+');
kIndex=dinp(1,1);
ik=dinp(2,1);

out_dir_name = strcat('results/case_k',num2str(kIndex),'_',num2str(ik));

lStyles='-';
colors=['g';'r';'m';'k'];
markers=['^';'o';'s';'v'];
lWidth=4;
mrkSize=18;

errWt=[1;1;1;1];

data=load('../data/SCW-683K.dat');
tData=data(1:5,1);
Data(:,1)=data(1:5,5)/100;
Data(:,2)=data(1:5,3)/100;
Data(:,3)=data(1:5,2)/100;
Data(:,4)=data(1:5,4)/100;
nTimes=numel(tData);

Mtot1=zeros(nTimes,4);

err=zeros(1,5);

errWtSum = ones(1,4)*errWt;
errWt = errWt/errWtSum;

fileName = strcat(out_dir_name,'/Time.dat');
dTime=load(fileName);
fileName = strcat(out_dir_name,'/Mr.dat');
Mr=load(fileName);
fileName = strcat(out_dir_name,'/Mt.dat');
Mt=load(fileName);
fileName = strcat(out_dir_name,'/Mck.dat');
dck=load(fileName);

t = dTime(:,1);
Mtot_tmp = Mr + Mt;
Mtot_tmp(:,end) = dck(:,3);    

Mtot(:,1)=Mtot_tmp(:,1);
Mtot(:,2)=Mtot_tmp(:,2) + Mtot_tmp(:,3) + Mtot_tmp(:,4) + Mtot_tmp(:,5) + Mtot_tmp(:,6);
Mtot(:,3)=Mtot_tmp(:,7) + Mtot_tmp(:,8);
Mtot(:,4)=Mtot_tmp(:,end);

for i=1:nTimes
    for j=1:numel(t)
        if(abs(t(j)-tData(i,1))<1e-3)
            for k=1:4
                Mtot1(i,k) = Mtot(j,k);
                err(1,k) = err(1,k) + (Mtot1(i,k) - Data(i,k))^2;
            end
            break;
        end
    end
end

err(1,5) = 0;
for k=1:4
    err(1,5) = err(1,5) + errWt(k,1)*err(1,k);
end

fig1 = figure(1);
hold on;
for k=1:4
    lStyleColor_tmp = strcat(lStyles,colors(k,1));
    markerColor_tmp = strcat(markers(k,1),colors(k,1));
    color_tmp = colors(k,1);

    plot(t(:,ik),Mtot(:,k,ik),lStyleColor_tmp,'linewidth',lWidth);
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

filename=sprintf('%s/err.dat',out_dir_name);
save(filename,'err','-ascii');