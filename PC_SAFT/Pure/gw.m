function spcsaftw1=gw(T,P,x,ro)
global m
global sigma
global epsk 
global kappa
global epsilon
global kij
nc=2;
R=8.314;
Na=6.02*10^23;
coef=(pi*Na)/6;
qe=[0 kij;kij 0];
tum0=0;
tum1=0;
tum2=0;
tum3=0;
mt=0;
for i=1:nc
  Bi(i)=epsk(i)/T;  
  d(i)=sigma(i)*(1-(0.12*exp(-3*Bi(i))));
  tum0=tum0+(x(i)*m(i));
  tum1=tum1+(x(i)*m(i)*d(i));
  tum2=tum2+(x(i)*m(i)*(d(i)^2));
  tum3=tum3+(x(i)*m(i)*(d(i)^3)); 
  mt=mt+(x(i)*m(i));
end
coef=(pi*Na)/6;
zeta0=ro*coef*tum0;
zeta1=ro*coef*tum1;
zeta2=ro*coef*tum2;
zeta3=ro*coef*tum3;
eta=ro*coef*tum3;
for i=1:nc
    for j=1:nc
        sigm(i,j)=0.5*(sigma(i)+sigma(j));
        epsdk(i,j)=(sqrt(epsk(i)*epsk(j)))*(1-qe(i,j));
    end
end
%--------------------------------------------------------------------------
Zhs=(6/(pi*Na*ro))*(((zeta0*zeta3)/(1-zeta3))+((3*zeta1*zeta2)/((1-zeta3)^2))+((3-zeta3)*(zeta2^3)/(1-zeta3)^3));
%--------------------------------------------------------------------------
Pum=0;
  for j=1:nc
      g(j)=(1/(1-zeta3))+((1.5*d(j))*(zeta2/(1-zeta3)^2))+(2*(0.5*d(j))^2)*((zeta2^2)/(1-zeta3)^3);
      dg(j)=(zeta3/(1-zeta3)^2)+(1.5*d(j)*zeta2/(1-zeta3)^2)+(3*d(j)*zeta2*zeta3/(1-zeta3)^3)+((d(j)^2)*(zeta2^2)/(1-zeta3)^3)+(1.5*(d(j)^2)*(zeta2^2)*zeta3/(1-zeta3)^4);
      Pum=Pum+(x(j)*(1-m(j))*(dg(j)/g(j)));
  end
Zchain=Pum;
%--------------------------------------------------------------------------
sum1=0;
sum2=0;
for i=1:nc
    for j=1:nc
      sum1=sum1+(m(i)*m(j)*x(i)*x(j)*(epsdk(i,j)/T)*(sigm(i,j)^3));
      sum2=sum2+(m(i)*m(j)*x(i)*x(j)*((epsdk(i,j)/T)^2)*(sigm(i,j)^3));
    end
end
con1=sum1;
con2=sum2;
%--------------------------------------------------------------------------
a0=[0.9105631445  0.6361281449   2.6861347891   -26.547362491  97.759208784   -159.59154087     91.297774084];    
a1=[-0.3084016918 0.1860531159   -2.5030047259   21.419793629  -65.255885330  83.318680481     -33.746922930];
a2=[-0.0906148351 0.4527842806   0.5962700728   -1.7241829131  -4.1302112531  13.776631870     -8.6728470368];
b0=[0.7240946941  2.2382791861   -4.0025849485  -21.003576815  26.855641363   206.55133841     -355.60235612];
b1=[-0.5755498075 0.6995095521   3.8925673390   -17.215471648  192.67226447   -161.82646165    -165.20769346];
b2=[0.0976883116  -0.2557574982  -9.1558561530   20.642075974  -38.804430052  93.626774077     -29.666905585];
%--------------------------------------------------------------------------
sum4=0;
sum5=0;
sum6=0;
sum7=0;
for i=1:7
    a(i)=a0(i)+(((mt-1)/mt)*a1(i))+((((mt-1)/mt)*((mt-2)/mt))*a2(i));
    b(i)=b0(i)+(((mt-1)/mt)*b1(i))+((((mt-1)/mt)*((mt-2)/mt))*b2(i));
    sum4=sum4+(a(i)*(eta^(i-1)));
    sum5=sum5+(b(i)*(eta^(i-1)));
    sum6=sum6+(i*a(i)*eta^(i-1));
    sum7=sum7+(i*b(i)*eta^(i-1));
end
I1=sum4;
I2=sum5;
dI1=sum6;
dI2=sum7;
%--------------------------------------------------------------------------
Par1=mt*(((8*eta)-(2*(eta^2)))/((1-eta)^4));
Par2=(1-mt)*((20*eta)-(27*(eta^2))+(12*(eta^3))-(2*(eta^4)));
Par3=((1-eta)*(2-eta))^2;
Para=1+Par1+(Par2/Par3);
%--------------------------------------------------------------------------
Tar1=4*mt*((2+(5*eta)-(eta^2))/((1-eta)^5));
Tar2=2*(1-mt)*((eta^3)+(6*(eta^2))-(24*eta)+20);
Tar3=((1-eta)*(2-eta))^3;
Tar=Tar1+(Tar2/Tar3);
Der=-Tar/(Para^2);
%--------------------------------------------------------------------------
s1=-pi*mt*ro*Na*(1/Para)*con2*dI2;
s2=-2*pi*ro*Na*con1*dI1;
s3=eta*Der*(-pi*mt*ro*Na*con2*I2);
Zdisp=s1+s2+s3;
%--------------------------------------------------------------------------
delta=(d(2)^3)*g(2)*kappa*(exp(epsilon/T)-1);
XA=(-1+sqrt(1+(4*Na*ro*x(2)*delta)))/(2*Na*ro*x(2)*delta);
ass1=1+(dg(2)/g(2));
Zassoc=-2*ass1*x(2)*(1-XA);
%--------------------------------------------------------------------------
spcsaftw1=(P/(ro*R*T))-1-Zhs-Zchain-Zdisp-Zassoc;
end
