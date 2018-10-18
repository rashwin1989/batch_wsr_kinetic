clear all;
close all;

nMrk=6;

dinp=load('plot_input.dat+');
kIndex=dinp(1,1);
ik=dinp(2,1);

out_dir_name = strcat('results/case_k',num2str(kIndex),'_',num2str(ik));

lStyles=['--';'--';'--';': ';'--';'--';'--';'--'];
colors=['g';'y';'y';'r';'r';'m';'m';'k'];
markers=['^';'>';'+';'o';'o';'+';'s';'v'];
lWidth=[4;4;4;4;4;4;4;4];
mrkSize=18;

% cmap1=colormap(autumn(5));
% orange=cmap1(2,:);
% cmap2=colormap(parula(5));
% teal=cmap2(3,:);
% cmap3=colormap(cool(5));
% dblue=cmap3(3,:);

out_dir_name = strcat('results/case_k',num2str(kIndex),'_',num2str(ik));

fileName=strcat(out_dir_name,'/Time.dat');
dTime=load(fileName);
fileName = strcat(out_dir_name,'/M1.dat');
Mo=load(fileName);
fileName = strcat(out_dir_name,'/M2.dat');
Mw=load(fileName);
fileName = strcat(out_dir_name,'/Mck.dat');
dck=load(fileName);

M1(:,1)=Mo(:,1);
M1(:,2)=Mo(:,2)+Mo(:,3);
M1(:,3)=Mo(:,4);
M1(:,4)=Mo(:,5);
M1(:,5)=Mo(:,6);
M1(:,6)=Mo(:,7);
M1(:,7)=Mo(:,8);
M1(:,8)=dck(:,1);

M2(:,1)=Mw(:,1);
M2(:,2)=Mw(:,2)+Mw(:,3);
M2(:,3)=Mw(:,4);
M2(:,4)=Mw(:,5);
M2(:,5)=Mw(:,6);
M2(:,6)=Mw(:,7);
M2(:,7)=Mw(:,8);
M2(:,8)=dck(:,2);

t = dTime(:,1);

fMrk = round(length(t)/nMrk);

fig1 = figure(1);    
hold on;
for k=1:numel(colors)
    lStyleColor_tmp = strcat(lStyles(k,:),colors(k,1));
    markerColor_tmp = strcat(markers(k,1),colors(k,1));
    color_tmp = colors(k,1);
    lWidth_tmp = lWidth(k,1);
    plot(t,M1(:,k),lStyleColor_tmp,'linewidth',lWidth_tmp);
    plot(t(1:fMrk:end),M1(1:fMrk:end,k),markerColor_tmp,'markersize',mrkSize, ...
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
filename=sprintf('%s/Mo.fig',out_dir_name);
savefig(fig1,filename);


fig2 = figure(2);    
hold on;
for k=1:numel(colors)
    lStyleColor_tmp = strcat(lStyles(k,:),colors(k,1));
    markerColor_tmp = strcat(markers(k,1),colors(k,1));
    color_tmp = colors(k,1);
    lWidth_tmp = lWidth(k,1);
    plot(t,M2(:,k),lStyleColor_tmp,'linewidth',lWidth_tmp);
    plot(t(1:fMrk:end),M2(1:fMrk:end,k),markerColor_tmp,'markersize',mrkSize, ...
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
filename=sprintf('%s/Mw.fig',out_dir_name);
savefig(fig2,filename);