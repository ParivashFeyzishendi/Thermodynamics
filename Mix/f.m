function bubble=f(coe)
global T
global m
global sigma
global epsk
global kappa
global epsilon
global nu
global xee
global kij
T=coe(1);
y(1)=coe(2);
y(2)=coe(3);
P=1.01325*100000;
nc=2;
%-------
%methanol && octanol
x(1)=0.0242;
x(2)=1-x(1);
Mw=[32.042 130.23];
%---------
R=8.314;
Na=6.02*10^23;
coef=(pi*Na)/6;
kij=0.020;
qe=[0 kij;kij 0];    
sigma=[3.2300*10^-10 3.7145*10^-10];
epsk=[188.90 262.74];
m=[1.5255 4.3555];
kappa=[0.035176 0.002197];
epsilon=[2899.5 2754.8];
%--------------------------------------------------------------------------
a0=[0.9105631445  0.6361281449   2.6861347891   -26.547362491  97.759208784   -159.59154087     91.297774084];    
a1=[-0.3084016918 0.1860531159   -2.5030047259   21.419793629  -65.255885330  83.318680481     -33.746922930];
a2=[-0.0906148351 0.4527842806   0.5962700728   -1.7241829131  -4.1302112531  13.776631870     -8.6728470368];
b0=[0.7240946941  2.2382791861   -4.0025849485  -21.003576815  26.855641363   206.55133841     -355.60235612];
b1=[-0.5755498075 0.6995095521   3.8925673390   -17.215471648  192.67226447   -161.82646165    -165.20769346];
b2=[0.0976883116  -0.2557574982  -9.1558561530   20.642075974  -38.804430052  93.626774077     -29.666905585];
%--------------------------------------------------------------------------
%Liquid
tuml0=0;
tuml1=0;
tuml2=0;
tuml3=0;
mlt=0;
for i=1:nc
  Bi(i)=epsk(i)/T;  
  d(i)=sigma(i)*(1-(0.12*exp(-3*Bi(i))));
  tuml0=tuml0+(x(i)*m(i));
  tuml1=tuml1+(x(i)*m(i)*d(i));
  tuml2=tuml2+(x(i)*m(i)*(d(i)^2));
  tuml3=tuml3+(x(i)*m(i)*(d(i)^3)); 
  mlt=mlt+(x(i)*m(i));
end
for i=1:nc
    for j=1:nc
        sigm(i,j)=0.5*(sigma(i)+sigma(j));
        epsdk(i,j)=(sqrt(epsk(i)*epsk(j)))*(1-qe(i,j));
    end
end
rol=6*(0.7)/(pi*Na*tuml3);
err=1;
while err>=0.001;
   megh=gw(T,P,x,rol);
   mosh=(gw(T,P,x,rol+(10^-5))-gw(T,P,x,rol-(10^-5)))/(2*(10^-5));
   ne=rol-(megh/mosh);
   err=abs(ne-rol);
   rol=ne;
end
molcw=0;
for i=1:nc
    molcw=molcw+(Mw(i)*x(i));
end
zetal0=rol*coef*tuml0;
zetal1=rol*coef*tuml1;
zetal2=rol*coef*tuml2;
zetal3=rol*coef*tuml3;
roLL=0.001*molcw*rol
etal=rol*coef*tuml3
ZL=P/(rol*R*T);
for i=1:nc
    rl(i)=rol*x(i);
end
%==========================================================================
%Liquid fugacity
%==========================================================================
% Hard sphere
con1=(3*zetal1*zetal2)/(1-zetal3)^2;
con2=zetal0/(1-zetal3);
con3=((zetal2^3)*((3*zetal3)-1))/((zetal3^2)*((1-zetal3)^3));
con4=log(1-zetal3);
con5=(2*con4)+(zetal3/(1-zetal3));
con6=con4+(zetal3/(1-zetal3)^2);
for i=1:nc
    Ql1(i)=(d(i)^3)*(con1+con2+con3);
    Ql2(i)=(3*d(i)*(zetal2+(d(i)*zetal1)))/(1-zetal3);
    Ql3(i)=(zetal2*d(i))/zetal3;
    Ql4(i)=(Ql3(i)^3)*con5;
    Ql5(i)=3*(Ql3(i)^2)*con6;
    muhsl(i)=m(i)*(Ql1(i)+Ql2(i)-con4-Ql4(i)+Ql5(i))
end
%==========================================================================
% Chain
for i=1:nc
      gl(i)=(1/(1-zetal3))+((1.5*d(i))*(zetal2/(1-zetal3)^2))+((2*(0.5*d(i))^2)*((zetal2^2)/(1-zetal3)^2));
end
for i=1:nc
   suml(i)=0;  
   for j=1:nc
     sl(j)=((d(i)^3)/(1-zetal3)^2)+((1.5*d(j)*(d(i)^2))/(1-zetal3)^2)+((3*d(j)*(d(i)^3)*zetal2)/(1-zetal3)^3)+(((d(j)^2)*(d(i)^2)*zetal2)/(1-zetal3)^3)+((1.5*(d(j)^2)*(d(i)^3)*(zetal2^2))/(1-zetal3)^4);      
     suml(i)=suml(i)+((rol*x(j))*(1-m(j))*(((pi/6)*m(i)*Na)*(sl(j)/gl(j))));
   end
   muchainl(i)=((1-m(i))*log(gl(i)))+suml(i)
end
%==========================================================================
% Dispersion
suml4=0;
suml5=0;
suml6=0;
suml7=0;
suml8=0;
suml9=0;
for i=1:7
    al(i)=a0(i)+(((mlt-1)/mlt)*a1(i))+((((mlt-1)/mlt)*((mlt-2)/mlt))*a2(i));
    bl(i)=b0(i)+(((mlt-1)/mlt)*b1(i))+((((mlt-1)/mlt)*((mlt-2)/mlt))*b2(i));
    suml4=suml4+(al(i)*(etal^(i-1)));
    suml5=suml5+(bl(i)*(etal^(i-1)));
    suml6=suml6+((((1/(mlt^2))*a1(i))+((((3*mlt)-4)/(mlt^3))*a2(i)))*(etal^(i-1)));
    suml7=suml7+((((1/(mlt^2))*b1(i))+((((3*mlt)-4)/(mlt^3))*b2(i)))*(etal^(i-1)));
    suml8=suml8+((i-1)*(etal^(i-2))*al(i));
    suml9=suml9+((i-1)*(etal^(i-2))*bl(i));
end
Il1=suml4;
Il2=suml5;
Dalm=suml6;
Dblm=suml7;
Deral=suml8;
Derbl=suml9;
suml1=0;
suml2=0;
for i=1:nc
    Addl1(i)=0;
    Addl2(i)=0;
    for j=1:nc
      suml1=suml1+(m(i)*m(j)*x(i)*x(j)*(epsdk(i,j)/T)*(sigm(i,j)^3));
      suml2=suml2+(m(i)*m(j)*x(i)*x(j)*((epsdk(i,j)/T)^2)*(sigm(i,j)^3));
      Addl1(i)=Addl1(i)+(2*m(i)*m(j)*x(j)*(epsdk(i,j)/T)*(sigm(i,j)^3));
      Addl2(i)=Addl2(i)+(2*m(i)*m(j)*x(j)*((epsdk(i,j)/T)^2)*(sigm(i,j)^3));
    end
end
conl1=suml1;
conl2=suml2;

for i=1:nc
    Detal(i)=coef*m(i)*rol*(d(i)^3);
    DIl1(i)=((m(i)-mlt)*Dalm)+(Detal(i)*Deral);
    DIl2(i)=((m(i)-mlt)*Dblm)+(Detal(i)*Derbl);
    Parl1(i)=(m(i)-mlt)*(((8*etal)-(2*(etal^2)))/((1-etal)^4));
    Parl2(i)=-(m(i)-mlt)*((20*etal)-(27*(etal^2))+(12*(etal^3))-(2*(etal^4)));
    Parl(i)=Parl1(i)+(Parl2(i)/(((1-etal)*(2-etal))^2));
end
Tarl1=4*mlt*((2+(5*etal)-(etal^2))/((1-etal)^5));
Tarl2=2*(1-mlt)*((etal^3)+(6*(etal^2))-(24*etal)+20);
Tarl3=((1-etal)*(2-etal))^3;
Tarl=Tarl1+(Tarl2/Tarl3);
Paral1=mlt*(((8*etal)-(2*(etal^2)))/((1-etal)^4));
Paral2=(1-mlt)*((20*etal)-(27*(etal^2))+(12*(etal^3))-(2*(etal^4)));
Paral3=((1-etal)*(2-etal))^2;
Paral=1+Paral1+(Paral2/Paral3);
for i=1:nc
  consl(i)=-((Detal(i)*Tarl)+Parl(i))/(Paral^2);
  Tsl1(i)=((m(i)-mlt)*(1/Paral)*Il2)+(mlt*consl(i)*Il2)+(mlt*(1/Paral)*DIl2(i));
  mudispl(i)=-(2*pi*rol*Na*Il1*Addl1(i))-(pi*rol*Na*(1/Paral)*mlt*Il2*Addl2(i))-(2*pi*rol*Na*DIl1(i)*conl1)-(pi*rol*Na*Tsl1(i)*conl2)
end
%--------------------------------------------------------------------------
xee=x;
nu=rol;
jav=fsolve(@fe,[0 0 0 0]);
XAl1=jav(1);
XBl1=jav(2);
XAl2=jav(3);
XBl2=jav(4);
for i=1:nc 
  for j=1:nc  
    sl(j,i)=((d(i)^3)/(1-zetal3)^2)+((1.5*d(j)*(d(i)^2))/(1-zetal3)^2)+((3*d(j)*(d(i)^3)*zetal2)/(1-zetal3)^3)+(((d(j)^2)*(d(i)^2)*zetal2)/(1-zetal3)^3)+((1.5*(d(j)^2)*(d(i)^3)*(zetal2^2))/(1-zetal3)^4);      
    assl(j,i)=(((pi/6)*m(i)*Na)*(sl(j,i)/gl(j)));
  end
 end
muassl(1)=(log(XAl1)+log(XBl1))-((0.5*rl(1)*((1-XAl1)+(1-XBl1))*assl(1,1))+(0.5*rl(2)*((1-XAl2)+(1-XBl2))*assl(2,1)))
muassl(2)=(log(XAl2)+log(XBl2))-((0.5*rl(1)*((1-XAl1)+(1-XBl1))*assl(1,2))+(0.5*rl(2)*((1-XAl2)+(1-XBl2))*assl(2,2)))
%--------------------------------------------------------------------------
%Vapor
tumv0=0;
tumv1=0;
tumv2=0;
tumv3=0;
mvt=0;
for i=1:nc
  tumv0=tumv0+(y(i)*m(i));
  tumv1=tumv1+(y(i)*m(i)*d(i));
  tumv2=tumv2+(y(i)*m(i)*(d(i)^2));
  tumv3=tumv3+(y(i)*m(i)*(d(i)^3)); 
  mvt=mvt+(y(i)*m(i));
end
rov=6*(10^-7)/(pi*Na*tuml3);
err=1;
while err>=0.001;
   megh=gw(T,P,y,rov);
   mosh=(gw(T,P,y,rov+(10^-5))-gw(T,P,y,rov-(10^-5)))/(2*(10^-5));
   ne=rov-(megh/mosh);
   err=abs(ne-rov);
   rov=ne;
end
molcv=0;
for i=1:nc
    molcv=molcv+(Mw(i)*y(i));
end
zetav0=rov*coef*tumv0;
zetav1=rov*coef*tumv1;
zetav2=rov*coef*tumv2;
zetav3=rov*coef*tumv3;
roVV=0.001*molcv*rov
etav=rov*coef*tumv3
ZV=P/(rov*R*T);
for i=1:nc
    rv(i)=rov*y(i);
end
%==========================================================================
%Vapor fugacity
%==========================================================================
% Hard sphere
cov1=(3*zetav1*zetav2)/(1-zetav3)^2;
cov2=zetav0/(1-zetav3);
cov3=((zetav2^3)*((3*zetav3)-1))/((zetav3^2)*((1-zetav3)^3));
cov4=log(1-zetav3);
cov5=(2*cov4)+(zetav3/(1-zetav3));
cov6=cov4+(zetav3/(1-zetav3)^2);
for i=1:nc
    Qv1(i)=(d(i)^3)*(cov1+cov2+cov3);
    Qv2(i)=(3*d(i)*(zetav2+(d(i)*zetav1)))/(1-zetav3);
    Qv3(i)=(zetav2*d(i))/zetav3;
    Qv4(i)=(Qv3(i)^3)*cov5;
    Qv5(i)=3*(Qv3(i)^2)*cov6;
    muhsv(i)=m(i)*(Qv1(i)+Qv2(i)-cov4-Qv4(i)+Qv5(i));
end
%==========================================================================
% Chain
for i=1:nc
      gv(i)=(1/(1-zetav3))+((1.5*d(i))*(zetav2/(1-zetav3)^2))+((2*(0.5*d(i))^2)*((zetav2^2)/(1-zetav3)^2));
end
for i=1:nc
   sumv(i)=0;  
   for j=1:nc
     sv(j)=((d(i)^3)/(1-zetav3)^2)+((1.5*d(j)*(d(i)^2))/(1-zetav3)^2)+((3*d(j)*(d(i)^3)*zetav2)/(1-zetav3)^3)+(((d(j)^2)*(d(i)^2)*zetav2)/(1-zetav3)^3)+((1.5*(d(j)^2)*(d(i)^3)*(zetav2^2))/(1-zetav3)^4);      
     sumv(i)=sumv(i)+((rov*y(j))*(1-m(j))*(((pi/6)*m(i)*Na)*(sv(j)/gv(j))));
   end
   muchainv(i)=((1-m(i))*log(gv(i)))+sumv(i);
end
%==========================================================================
% Dispersion
sumv4=0;
sumv5=0;
sumv6=0;
sumv7=0;
sumv8=0;
sumv9=0;
for i=1:7
    av(i)=a0(i)+(((mvt-1)/mvt)*a1(i))+((((mvt-1)/mvt)*((mvt-2)/mvt))*a2(i));
    bv(i)=b0(i)+(((mvt-1)/mvt)*b1(i))+((((mvt-1)/mvt)*((mvt-2)/mvt))*b2(i));
    sumv4=sumv4+(av(i)*(etav^(i-1)));
    sumv5=sumv5+(bv(i)*(etav^(i-1)));
    sumv6=sumv6+((((1/(mvt^2))*a1(i))+((((3*mvt)-4)/(mvt^3))*a2(i)))*(etav^(i-1)));
    sumv7=sumv7+((((1/(mvt^2))*b1(i))+((((3*mvt)-4)/(mvt^3))*b2(i)))*(etav^(i-1)));
    sumv8=sumv8+((i-1)*(etav^(i-2))*av(i));
    sumv9=sumv9+((i-1)*(etav^(i-2))*bv(i));
end
Iv1=sumv4;
Iv2=sumv5;
Davm=sumv6;
Dbvm=sumv7;
Derav=sumv8;
Derbv=sumv9;
sumv1=0;
sumv2=0;
for i=1:nc
    Addv1(i)=0;
    Addv2(i)=0;
    for j=1:nc
      sumv1=sumv1+(m(i)*m(j)*y(i)*y(j)*(epsdk(i,j)/T)*(sigm(i,j)^3));
      sumv2=sumv2+(m(i)*m(j)*y(i)*y(j)*((epsdk(i,j)/T)^2)*(sigm(i,j)^3));
      Addv1(i)=Addv1(i)+(2*m(i)*m(j)*y(j)*(epsdk(i,j)/T)*(sigm(i,j)^3));
      Addv2(i)=Addv2(i)+(2*m(i)*m(j)*y(j)*((epsdk(i,j)/T)^2)*(sigm(i,j)^3));
    end
end
conv1=sumv1;
conv2=sumv2;

for i=1:nc
    Detav(i)=coef*m(i)*rov*(d(i)^3);
    DIv1(i)=((m(i)-mvt)*Davm)+(Detav(i)*Derav);
    DIv2(i)=((m(i)-mvt)*Dbvm)+(Detav(i)*Derbv);
    Parv1(i)=(m(i)-mvt)*(((8*etav)-(2*(etav^2)))/((1-etav)^4));
    Parv2(i)=-(m(i)-mvt)*((20*etav)-(27*(etav^2))+(12*(etav^3))-(2*(etav^4)));
    Parv(i)=Parv1(i)+(Parv2(i)/(((1-etav)*(2-etav))^2));
end
Tarv1=4*mvt*((2+(5*etav)-(etav^2))/((1-etav)^5));
Tarv2=2*(1-mvt)*((etav^3)+(6*(etav^2))-(24*etav)+20);
Tarv3=((1-etav)*(2-etav))^3;
Tarv=Tarv1+(Tarv2/Tarv3);
Parav1=mvt*(((8*etav)-(2*(etav^2)))/((1-etav)^4));
Parav2=(1-mvt)*((20*etav)-(27*(etav^2))+(12*(etav^3))-(2*(etav^4)));
Parav3=((1-etav)*(2-etav))^2;
Parav=1+Parav1+(Parav2/Parav3);
for i=1:nc
   consv(i)=-((Detav(i)*Tarv)+Parv(i))/(Parav^2);
   Tsv1(i)=((m(i)-mvt)*(1/Parav)*Iv2)+(mvt*consv(i)*Iv2)+(mvt*(1/Parav)*DIv2(i));
   mudispv(i)=-(2*pi*rov*Na*Iv1*Addv1(i))-(pi*rov*Na*(1/Parav)*mvt*Iv2*Addv2(i))-(2*pi*rov*Na*DIv1(i)*conv1)-(pi*rov*Na*Tsv1(i)*conv2);
end
%--------------------------------------------------------------------------
xee=y;
nu=rov;
javv=fsolve(@fe,[0 0 0 0]);
XAv1=javv(1);
XBv1=javv(2);
XAv2=javv(3);
XBv2=javv(4);
for i=1:nc 
  for j=1:nc  
    sv(j,i)=((d(i)^3)/(1-zetav3)^2)+((1.5*d(j)*(d(i)^2))/(1-zetav3)^2)+((3*d(j)*(d(i)^3)*zetav2)/(1-zetav3)^3)+(((d(j)^2)*(d(i)^2)*zetav2)/(1-zetav3)^3)+((1.5*(d(j)^2)*(d(i)^3)*(zetav2^2))/(1-zetav3)^4);      
    assv(j,i)=(((pi/6)*m(i)*Na)*(sv(j,i)/gv(j)));
  end
 end
muassv(1)=(log(XAv1)+log(XBv1))-((0.5*rv(1)*((1-XAv1)+(1-XBv1))*assv(1,1))+(0.5*rv(2)*((1-XAv2)+(1-XBv2))*assv(2,1)))
muassv(2)=(log(XAv2)+log(XBv2))-((0.5*rv(1)*((1-XAv1)+(1-XBv1))*assv(1,2))+(0.5*rv(2)*((1-XAv2)+(1-XBv2))*assv(2,2)))
%--------------------------------------------------------------------------
for i=1:nc
  muresl(i)=muhsl(i)+muchainl(i)+mudispl(i)+muassl(i)
  muresv(i)=muhsv(i)+muchainv(i)+mudispv(i)+muassv(i)
  phil(i)=exp(muresl(i)-log(ZL));
  phiv(i)=exp(muresv(i)-log(ZV))
  fl(i)=P*phil(i)*x(i)
  fv(i)=P*phiv(i)*y(i)
end
bubble=[(fl(1)/fv(1))-1
        (fl(2)/fv(2))-1
        y(1)+y(2)-1]
end
