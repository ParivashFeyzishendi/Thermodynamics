format long g
clear
clc
T=100;
Tc=190.4;
Pc=46*101325;
w=0.0115;
Tr=T/Tc;
f00=5.92714-(6.09648/Tr)-1.28862*log(Tr)+0.16934*(Tr)^6;
f1=15.2518-(15.6875/Tr)-13.4721*log(Tr)+0.43577*(Tr)^6;
P=Pc*(exp(f00+(w*f1)));
sigma=3.7039*10^-10;
epsk=150.03;
m=1.0000;
Bi=epsk/T;
Na=6.02*10^23;
R=8.314;
d=sigma*(1-(0.12*exp(-3*Bi)));
Eror=4;
while (Eror>=0.001)
rol=(6*(0.7))/(pi*Na*m*(d^3));
Er=1;
while Er>10^-8
   megh=f(T,P,rol,m,sigma,epsk);
   mosh=(f(T,P,rol+(10^-5),m,sigma,epsk)-f(T,P,rol-(10^-5),m,sigma,epsk))/(2*(10^-5));
   ne=rol-(megh/mosh);
   Er=abs(ne-rol);
   rol=ne;  
end
rol
Zl=P/(rol*R*T);
etal=((pi*Na)/6)*rol*m*(d^3)
rov=(6*(1e-7))/(pi*Na*m*(d^3));
Er=1;
while Er>10^-8
   megh=f(T,P,rov,m,sigma,epsk);
   mosh=(f(T,P,rov+(10^-5),m,sigma,epsk)-f(T,P,rov-(10^-5),m,sigma,epsk))/(2*(10^-5));
   ne=rov-(megh/mosh);
   Er=abs(ne-rov);
   rov=ne;  
end
rov
Zv=P/(rov*R*T);
etav=((pi*Na)/6)*rov*m*(d^3);
%--------------------------------------------------------------------------
a0=[0.9105631445  0.6361281449   2.6861347891   -26.547362491  97.759208784   -159.59154087     91.297774084];    
a1=[-0.3084016918 0.1860531159   -2.5030047259   21.419793629  -65.255885330  83.318680481     -33.746922930];
a2=[-0.0906148351 0.4527842806   0.5962700728   -1.7241829131  -4.1302112531  13.776631870     -8.6728470368];
b0=[0.7240946941  2.2382791861   -4.0025849485  -21.003576815  26.855641363   206.55133841     -355.60235612];
b1=[-0.5755498075 0.6995095521   3.8925673390   -17.215471648  192.67226447   -161.82646165    -165.20769346];
b2=[0.0976883116  -0.2557574982  -9.1558561530   20.642075974  -38.804430052  93.626774077     -29.666905585];
%--------------------------------------------------------------------------
Ahsv=R*T*m*(((4*etav)-(3*(etav^2)))/((1-etav)^2));
Achainv=R*T*(1-m)*log((1-(0.5*etav))/((1-etav)^3));
Parv1=m*(((8*etav)-(2*(etav^2)))/((1-etav)^4));
Parv2=(1-m)*((20*etav)-(27*(etav^2))+(12*(etav^3))-(2*(etav^4)));
Parv3=((1-etav)*(2-etav))^2;
Parav=1+Parv1+(Parv2/Parv3);
sumv1=0;
sumv2=0;
for i=1:7
    a(i)=a0(i)+(((m-1)/m)*a1(i))+((((m-1)/m)*((m-2)/m))*a2(i));
    b(i)=b0(i)+(((m-1)/m)*b1(i))+((((m-1)/m)*((m-2)/m))*b2(i));
    sumv1=sumv1+(a(i)*(etav^(i-1)));
    sumv2=sumv2+(b(i)*(etav^(i-1)));
end
Iv1=sumv1;
Iv2=sumv2;
A1v=-2*pi*rov*Na*Bi*(m^2)*(sigma^3)*Iv1;
A2v=-pi*rov*m*Na*(1/Parav)*(Bi^2)*(m^2)*(sigma^3)*Iv2;
Adispv=R*T*(A1v+A2v);
%--------------------------------------------------------------------------
Ahsl=R*T*m*(((4*etal)-(3*(etal^2)))/((1-etal)^2));
Achainl=R*(1-m)*T*log((1-(0.5*etal))/((1-etal)^3));
Parl1=m*(((8*etal)-(2*(etal^2)))/((1-etal)^4));
Parl2=(1-m)*((20*etal)-(27*(etal^2))+(12*(etal^3))-(2*(etal^4)));
Parl3=((1-etal)*(2-etal))^2;
Paral=1+Parl1+(Parl2/Parl3);
suml1=0;
suml2=0;
for i=1:7
    a(i)=a0(i)+(((m-1)/m)*a1(i))+((((m-1)/m)*((m-2)/m))*a2(i));
    b(i)=b0(i)+(((m-1)/m)*b1(i))+((((m-1)/m)*((m-2)/m))*b2(i));
    suml1=suml1+(a(i)*(etal^(i-1)));
    suml2=suml2+(b(i)*(etal^(i-1)));
end
Il1=suml1;
Il2=suml2;
A1l=-2*pi*rol*Na*Bi*(m^2)*(sigma^3)*Il1;
A2l=-pi*rol*m*Na*(1/Paral)*(Bi^2)*(m^2)*(sigma^3)*Il2;
Adispl=R*T*(A1l+A2l);
%--------------------------------------------------------------------------
Aresv=Ahsv+Achainv+Adispv;
Aresl=Ahsl+Achainl+Adispl;
SHv=(Aresv/(R*T))+(Zv-1)-log(Zv);
Phiv=exp(SHv);
SHl=(Aresl/(R*T))+(Zl-1)-log(Zl);
Phil=exp(SHl);
Pnew=P*(Phil/Phiv);
Eror=abs(Pnew-P);
P=Pnew;
end
Peq=P/1013250
rol=rol




