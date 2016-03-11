function range_4
D=2.4E-9;  %m^2/s
Na=25E-6;      %mol/m^3-s, if positive, get a positive flux,
               %negative concentration values, opposite for negative
%vary N, keep D = 2.4E-9 as stated earlier
C1 = mass_transfer(0, Na);
C2 = mass_transfer(D, Na);
C3 = mass_transfer(D*10, Na);
C4 = mass_transfer(D*100, Na);
plot_graphs(C1, C2, C3, C4, D, Na)


function Conc = mass_transfer(D, Na)
 
close all

Nx=201;    %number of points in x direction
Nt=101;    %number of points in time
m=0;       %cartesian
%D=2.4E-9;  %m^2/s
%Na=-25E-6;      %mol/m^3-s, if positive, get a positive flux,
               %negative concentration values, opposite for negative
Cw=4E-8;   %mol/m^3
Ci=4E-8;      %mol/m^3
w=40E-6;    %m
 
Xinf=.1;   %m
tmax=200000;  %sec
dx=Xinf/(Nx-1);
dt=tmax/(Nt-1);
 
t=linspace(0,tmax,Nt);
x=linspace(0,Xinf,Nx);
options=odeset;  %('MaxStep',100);
pbBC=@(yl,Cl,yr,Cr,t)pbBC1(yl,Cl,yr,Cr,t,Cw);
pdepb=@(x,t,C,DCDx)pdepb1(x,t,C,DCDx,Na,D);
sol=pdepe(m,pdepb,@pbIC,pbBC,x,t,options);  %solution      
C=sol(:,:,1);                          %extract solution 
Conc=C;
 
function res = plot_graphs(C1, C2, C3, C4, D, Na)
Nx=201;    %number of points in x direction
Nt=101;    %number of points in time
m=0;       %cartesian
%D=2.4E-9;  %m^2/s
%Na=-25E-6;      %mol/m^3-s, if positive, get a positive flux,
               %negative concentration values, opposite for negative
Cw=4E-8;   %mol/m^3
Ci=0;      %mol/m^3
w=40E-6;    %m
 
Xinf=.1;   %m
tmax=200000;  %sec
dx=Xinf/(Nx-1);
dt=tmax/(Nt-1);

t=linspace(0,tmax,Nt);
x=linspace(0,Xinf,Nx);

x_value1=6;
x_value2=16;   
x_value3=36; 
x_value4=51;
plot(t,C1(:,x_value1),t,C1(:,x_value2),t,C1(:,x_value3),t,C1(:,x_value4)) %Plot vs t at specified x locations
hold on 
plot(t,C2(:,x_value1),':',t,C2(:,x_value2),':',t,C2(:,x_value3),':',t,C2(:,x_value4),':')   %Plot vs t at specified x locations
plot(t,C3(:,x_value1),'--',t,C3(:,x_value2),'--',t,C3(:,x_value3),'--',t,C3(:,x_value4),'--')   %Plot vs t at specified x locations
plot(t,C4(:,x_value1),'-.',t,C4(:,x_value2),'-.',t,C4(:,x_value3),'-.',t,C4(:,x_value4),'-.')   %Plot vs t at specified x locations
hold off
xlabel('Time (s)')
ylabel('Concentration (mol/m^3)')
legend([num2str(x_value1*dx-dx) ' m'],[num2str(x_value2*dx-dx) ' m'],[num2str(x_value3*dx-dx) ' m'],[num2str(x_value4*dx-dx) ' m'],['- = 0 mol/m^3-s'],['-- = 2.4E-9 mol/m^3-s'],[': = 2.4E-10 mol/m^3-s'], ['-. = 2.4E-11 mol/m^3-s'])
 
figure(2)   %Plot vs x at specified times
t1=2;
t2=11;
t3=31;
t4=101;
plot(x,C1(t1,:),x,C1(t2,:),x,C1(t3,:),x,C1(t4,:))
hold on 
plot(x,C2(t1,:),':',x,C2(t2,:),':',x,C2(t3,:),':',x,C2(t4,:),':')
plot(x,C3(t1,:),'--',x,C3(t2,:),'--',x,C3(t3,:),'--',x,C3(t4,:),'--')
plot(x,C4(t1,:),'-.',x,C4(t2,:),'-.',x,C4(t3,:),'-.',x,C4(t4,:),'-.')
hold off
xlabel('x (m)')
ylabel('Concentration (mol/m^3)')
%ylim([0 Cw*(40)])
xlim([0 w*(1E-3)])
legend([num2str(t1*dt-dt) ' sec'],[num2str(t2*dt-dt) ' sec'],[num2str(t3*dt-dt) ' sec'], [num2str(t4*dt-dt) ' sec'])
 
%flux at base (x=0)
figure(3)
dCdx=(C1(:,2)-C1(:,1))/dx;
flux1=-D*dCdx;
plot(t,flux1)
hold on 
dCdx=(C2(:,2)-C2(:,1))/dx;
flux2=-D*dCdx;
plot(t,flux2,'--')
dCdx=(C3(:,2)-C3(:,1))/dx;
flux3=-D*dCdx;
plot(t,flux3,':')
dCdx=(C4(:,2)-C4(:,1))/dx;
flux4=-D*dCdx;
plot(t,flux4,'-.')
hold off
xlabel('Time (sec)')
ylabel('Flux (mol/s-m^2)')
figure(4)
loglog(t,flux1)
hold on 
loglog(t,flux2,'--')
loglog(t,flux3,':')
loglog(t,flux4,'-.')
hold off
xlabel('Time (sec)')
ylabel('Flux (mol/s-m^2)')


function [c,f,s]=pdepb1(x,t,C,DCDx,Na,D)  %define pde system – see Matlab help for pdepe
c=1;           %coefficient on time derivative
f=D*DCDx;      %flux   mol/s-m^3
s=-Na;         %source   mol/s-m^3
 
function C0=pbIC(x)   %initial condition
C0=4E-8;   %mol/m^3
 
function [pl,ql,pr,qr]=pbBC1(yl,Cl,yr,Cr,t,Cw)  %boundary conditions
pl=Cl;   %left side bc. This one sets the concentration to 0 (impermeable)
ql=0;       %Form is pl + ql*f = 0 where f is flux from definition of 
            %pde system in function pdepb
pr=Cr-Cw;      %right side bc. This one sets the concentration to Cw
qr=0;