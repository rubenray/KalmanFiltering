load 'datosKal.mat' 
%%time stamp,Temp Aire, Temp interna; Radiacion, Tem Aire 2 la ultima solo para referencia
% carga prueba 1 , pueba 2 , caraga temp,carga2,kalmanfusion2,loadValdez
%% dT/dt= P/mc - Kr(Tin-Ta); Kr=ha/mc P=potencia entrada m masa c water specific heat 
%%Kr  Kr=- dT/(dt*(T-Ta)) con P=0;Kr 
data=datosKal(1:4,:);
dT=(data(1,2)-data(1,1)); %T=30 secs  3.472222015261650e-04 needed to compute power
T=30;
data(2,:)=datosKal(5,:);
%time Tem_amb, Tem_int, Rad '
%uncomment below lines to compute sampling rate in
% this case 30 seg
% datestr(data(1,1)) ; 
% datestr(data(1,2))

[numobs,LL]=size(data');

numstates=6;

for k=2:1:40000
    dderivative(k)=(2*data(3,k)-data(3,k-1)-data(3,k+1))/(T^2);
end
plot(dderivative(100:1000))

