function comp=f(T,P,ro,m,sigma,epsk)
R=8.314;
Na=6.02*10^23;
Bi=epsk/T;
d=sigma*(1-(0.12*exp(-3*Bi)));
eta=((pi*Na)/6)*ro*m*(d^3);
Zhs=m*(((4*eta)-(2*(eta^2)))/((1-eta)^3));
Zchain=(1-m)*(((5*eta)-(2*(eta^2)))/((1-eta)*(2-eta)));
%--------------------------------------------------------------------------
Par1=m*(((8*eta)-(2*(eta^2)))/((1-eta)^4));
Par2=(1-m)*((20*eta)-(27*(eta^2))+(12*(eta^3))-(2*(eta^4)));
Par3=((1-eta)*(2-eta))^2;
Para=1+Par1+(Par2/Par3);
%--------------------------------------------------------------------------
Tar1=4*m*((2+(5*eta)-(eta^2))/((1-eta)^5));
Tar2=2*(1-m)*((eta^3)+(6*(eta^2))-(24*eta)+20);
Tar3=((1-eta)*(2-eta))^3;
Tar=Tar1+(Tar2/Tar3);
Der=-Tar/(Para^2);
%--------------------------------------------------------------------------
a0=[0.9105631445  0.6361281449   2.6861347891   -26.547362491  97.759208784   -159.59154087     91.297774084];    
a1=[-0.3084016918 0.1860531159   -2.5030047259   21.419793629  -65.255885330  83.318680481     -33.746922930];
a2=[-0.0906148351 0.4527842806   0.5962700728   -1.7241829131  -4.1302112531  13.776631870     -8.6728470368];
b0=[0.7240946941  2.2382791861   -4.0025849485  -21.003576815  26.855641363   206.55133841     -355.60235612];
b1=[-0.5755498075 0.6995095521   3.8925673390   -17.215471648  192.67226447   -161.82646165    -165.20769346];
b2=[0.0976883116  -0.2557574982  -9.1558561530   20.642075974  -38.804430052  93.626774077     -29.666905585];
%--------------------------------------------------------------------------
sum1=0;
sum2=0;
for i=1:7
    a(i)=a0(i)+(((m-1)/m)*a1(i))+((((m-1)/m)*((m-2)/m))*a2(i));
    b(i)=b0(i)+(((m-1)/m)*b1(i))+((((m-1)/m)*((m-2)/m))*b2(i));
    sum1=sum1+(a(i)*(eta^(i-1)));
    sum2=sum2+(b(i)*(eta^(i-1)));
end
I1=sum1;
I2=sum2;
%--------------------------------------------------------------------------
sum3=0;
sum4=0;
for i=1:7
  sum3=sum3+(i*a(i)*eta^(i-1));
  sum4=sum4+(i*b(i)*eta^(i-1));
end
%--------------------------------------------------------------------------
dI1=sum3;
dI2=sum4;
%--------------------------------------------------------------------------
s1=-pi*m*ro*Na*(1/Para)*(Bi^2)*(m^2)*(sigma^3)*dI2;
s2=-2*pi*ro*Na*Bi*(m^2)*(sigma^3)*dI1;
s3=eta*Der*(-pi*m*ro*Na*(Bi^2)*(m^2)*(sigma^3)*I2);
Zdisp=s1+s2+s3;
%--------------------------------------------------------------------------
comp=(P/(ro*R*T))-1-Zhs-Zchain-Zdisp;
end
