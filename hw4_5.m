function OneD_transient_mass_transfer
 
close all

 
Nx=201;    %number of points in x direction
Nt=101;    %number of points in time
m=0;       %cartesian
D=2.4E-9;  %m^2/s
Na=0;      %mol/m^3-s, if positive, get a positive flux,
               %negative concentration values, opposite for negative
Cw=4E-8;   %mol/m^3
%Ci=0;      %mol/m^3
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

x_value1=6;
x_value2=16;   
x_value3=36; 
x_value4=51;
plot(t,C(:,x_value1),t,C(:,x_value2),t,C(:,x_value3),t,C(:,x_value4))   %Plot vs t at specified x locations
xlabel('Time (s)')
ylabel('Concentration (mol/m^3)')
legend([num2str(x_value1*dx-dx) ' m'],[num2str(x_value2*dx-dx) ' m'],[num2str(x_value3*dx-dx) ' m'],[num2str(x_value4*dx-dx) ' m'])
 
figure(2)   %Plot vs x at specified times
t1=2;
t2=11;
t3=31;
t4=101;
plot(x,C(t1,:),x,C(t2,:),x,C(t3,:),x,C(t4,:))
xlabel('x (m)')
ylabel('Concentration (mol/m^3)')
%ylim([0 Cw*(40)])
%xlim([0 w*(1E-3)])
legend([num2str(t1*dt-dt) ' sec'],[num2str(t2*dt-dt) ' sec'],[num2str(t3*dt-dt) ' sec'], [num2str(t4*dt-dt) ' sec'])
 
%flux at base (x=0)
figure(3)
dCdx=(C(:,2)-C(:,1))/dx;
flux=-D*dCdx;
plot(t,flux)
xlabel('Time (sec)')
ylabel('Flux (mol/s-m^2)')
figure(4)
loglog(t,flux)
xlabel('Time (sec)')
ylabel('Flux (mol/s-m^2)')


function [c,f,s]=pdepb1(x,t,C,DCDx,Na,D)  %define pde system – see Matlab help for pdepe
c=1;           %coefficient on time derivative
f=D*DCDx;      %flux   mol/s-m^3
s=-Na;         %source   mol/s-m^3
 
function C0=pbIC(x)   %initial condition, this parameter was varied 
                      %throughout each run 
C0=sin(x/2);   %mol/m^3
 
function [pl,ql,pr,qr]=pbBC1(yl,Cl,yr,Cr,t,Cw)  %boundary conditions
pl=Cl-Cw;   %left side bc. This one sets the concentration to 0 (impermeable)
ql=0;       %Form is pl + ql*f = 0 where f is flux from definition of 
            %pde system in function pdepb
pr=Cr;      %right side bc. This one sets the concentration to Cw
qr=0;