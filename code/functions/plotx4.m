function [hh,hh1,hh2]=plotx4(t,y)
set(gcf,'DefaultAxesColorOrder',[0.8 0.1 0.1;1 0 0;1 0 0;0 0 1]);
cu=y(:,2);
cl=y(:,3);
du=y(:,4);
dl=y(:,5);

h=t;
h=h;
hh=fill([h(1); h(1:end); flipud([h(1:end); h(end)])],[cu(1); cl(1:end); flipud([cu(1:end); cl(size(cl,1))])],'b');
set(hh,'edgecolor',[1 0.75 0.75]);
set(hh,'facecolor',[1 0.75 0.75]);

hold on;

hh=fill([h(1); h(1:end); flipud([h(1:end); h(end)])],[du(1); dl(1:end); flipud([du(1:end); dl(size(cl,1))])],'b');
set(hh,'edgecolor',[1 0.90 0.90]);
set(hh,'facecolor',[1 0.90 0.90]);
hold on
 hh1=plot(h,y(:,1),'r','LineWidth',2);

hold on;
zz=zeros(size(y,1),1);
hh2=plot(h,zz,'b-');
