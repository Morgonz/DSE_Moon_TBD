n=100;
B=zeros(n,3);
for i=(1:n)
    C=MPPVfun(10000+1737+1629-10000/n*i);
    B(i,:)=C;
end

ax1=subplot(2,2,1);
ax2=subplot(2,2,2);
ax3=subplot(2,2,3);
p1=plot(ax1,B(:,1),B(:,2));
p2=plot(ax2,B(:,1),B(:,3));
p3=plot(ax3,B(:,2),B(:,3));
%axis(ax1,[1629+1737,13000]);
grid(ax1,'on');
grid(ax1,'minor');
%axis(ax2,[1629+1737,13000]);
grid(ax2,'on');
grid(ax2,'minor');
grid(ax3,'on');
grid(ax3,'minor');
title(ax1, '2 dep. vehicles')
xlabel(ax1,'apoapse[km]');
xlabel(ax2,'apoapse[km]');
xlabel(ax3,'Dep. mass [km]');
ylabel(ax3,'Sat. mass [km]');
ylabel(ax1,'Dep. mass [kg]');
ylabel(ax2,'Sat. mass [kg]');
%p1.Marker='x';
%p2.Marker='x';