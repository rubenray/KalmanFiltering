load 'datosKal.mat'
%time stamp,Temp Aire, Temp interna; Radiacion, Tem Aire2
%
%% dT/dt= P/mc - Kr(Tin-Ta); Kr=ha/mc P=power m mass c water specific heat
%% Kr=- dT/(dt*(T-Ta)) con P=0;Kr Transfer coheficient
data=datosKal(1:4,:);
dT=(data(1,2)-data(1,1)); %T=30 secs   needed to compute power
T=30;
data(2,:)=datosKal(5,:);
%uncomment below lines to compute sampling rate in
% this case 30 seg
% datestr(data(1,1)) ;
% datestr(data(1,2))

[numobs,L]=size(data');
numstates=6;
p11=zeros(numobs,1);
p22=zeros(numobs,1);
p33=zeros(numobs,1);
p44=zeros(numobs,1);
p55=zeros(numobs,1);
p66=zeros(numobs,1);
qf=.001;
qs=.000001;%qs=.001
qp1=1e-10;
qp2=1e-14;
Kr=(29.19-28.44)/[(7*3600)*(29.18-23)];

%A=pi*r^2
%A=pi*(1.55/2)^2;
A=1.9;%m^2

c =4181; %joules/kgm
volumen =(1.75)*(1.9); %3.32^3
%mass=rho*vol ~~~ 1000 ;
m=1000;
% u1=zeros(numobs,1);
% u2=zeros(numobs,1);
% u3=zeros(numobs,1);
u4=data(4,:)';
u=u4;
Kra=Kr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=[1 0 0 0 0 0; 0 0 1 0 0  0 ] ;

P=[ 100 0  0  0  0  0; 0  100 0  0  0  0; 0  0  10 0 0 0; 0  0   0 10  0  0; 0 0 0 0 10 0; 0 0 0 0 0 .1];
r1=0.01;
r2=1;

R=[r1 0; 0 r2];
Q1=[T^3*qf  T^2*qf ; T^2*qf    T*qf  ];
Q2=[0     0; 0     0];
Q3=[0     0; 0      0];
Q4=[T^3*qs  T^2*qs ; T^2*qs    T*qs  ];
Qz=zeros(2,2);
Qp=[qp1  0 ; 0  qp2];

Q= [T^3*qf  T^2*qf       0       0      0      1e-10 ;
    T^2*qf    T*qf       0       0      0       0  ;
    0           0     T^3*qs  T^2*qs    1e-10  1e-9;
    0           0     T^2*qs    T*qs    0       0  ;
    0           0    1e-10      0       qp1       0  ;
    1e-10      0    1e-9       0       0      qp2];


% G=5.4532e-08;
% G=5.8e-08
% G=.0016
% G=6.47e-05
% G=2.67e-6
% G=1.63e-07
G=6.24e-08;
%G=5.8e-08

Kr=Kra/30;

%numobs=40000;
xhat=zeros(numstates,numobs);
xhat(:,1)=[data(2,1),0,data(3,1), 0,Kr ,G ];
xh=xhat(:,1);

%%%%kalman process
for k = 2:numobs
    
    Kr=xhat(5,k-1);
    G=xhat(6,k-1);
    
    F=   [1    T    0     0     0   0  ;
        0    1    0     0     0   0;
        0    0    1     T     0   0;
        +Kr*T 0  -Kr*T   0     0   0 ;
        0   0    0      0     1   0;
        0   0    0      0     0   1;
        ];
    
    B=[0 ,0 ,0, G , 0 ,  0]';
    
    xh= F*xh+B*u(k);%0
    
    xh(1,1)=data(2,k) ;
    
    P=F*P*F'+Q;
    
    K    = P  * H' / (H * P * H' + R);
    
    P    = P - K * H * P;
    p11(k)=P(1,1); p22(k)=P(2,2);p33=P(3,3);p44=P(4,4);p55=P(5,5);p66=P(6,6);
    
    s1=data(2,k); s2=data(3,k);
    
    if isnan(s1)
        s1=xh(1,1);
    end
    if isnan(s2)
        s2=xh(3,1);
    end
    k;
    z=[s1; s2];
    xh= xh+K*(z -H*xh);
    xhat(:,k)=xh;
    %xhat(5,k) uncomment to plot xhat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumTicks = 12;
figure(1);
subplot(2,2,1);
plot(data(1,1:numobs),data(3,1:numobs),'-k');
title('Internal Temperature with missing points ')
grid on
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
datetick('x','mm/dd','keeplimits', 'keepticks')
xticklabel_rotate;

subplot(2,2,2);
plot(data(1,1:numobs),xhat(3,1:numobs),'-b','LineWidth',2);
title('Filtered Internal Temperature')
grid on
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
datetick('x','mm/dd','keeplimits', 'keepticks')
xticklabel_rotate;

subplot(2,2,3);
plot(data(1,1:numobs),data(3,1:numobs),'or');
hold on
plot(data(1,1:numobs),xhat(3,1:numobs),'-b','LineWidth',2);

%plot(data(:,1),xhatT(2,:)');
%plot(data(:,1),xhatT(4,:)');
%plot(data(1,1:numobs),xhat(1,1:numobs));
%plot(data(1,1:numobs),u(1:numobs)/100+20);

title('Original plus Filtered Signal ')
grid on
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
datetick('x','mm/dd','keeplimits', 'keepticks')
xticklabel_rotate;


subplot(2,2,4);
%plot(data(1,:),xhat(3,1:numobs),'or');
plot(data(1,:),xhat(1,:),'LineWidth',2);
hold on
plot(data(1,:),u/50,'LineWidth',2);
plot(data(1,:),xhat(3,:),'LineWidth',2 );
title('External and Internal Temperature and Radition')
grid on
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
datetick('x','mm/dd','keeplimits', 'keepticks')
xticklabel_rotate;
hold off;

figure(2)
plot(xhat(5,:));
set(gca,'FontSize',14);
title('1  Parameter estimation');
grid on

figure(3)
plot(xhat(6,:));
title('2  Parameter estimation')
set(gca,'FontSize',14);
grid on
