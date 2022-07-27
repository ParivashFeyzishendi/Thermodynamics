function moadele=fe(coe)
global T
global m
global sigma
global epsk
global kappa
global epsilon
global nu
global xee
x(1)=xee(1);
x(2)=xee(2);
Na=6.02*10^23;
coef=(pi*Na)/6;
n=2;
tum2=0;
tum3=0;
for i=1:n
  Bi(i)=epsk(i)/T;  
  d(i)=sigma(i)*(1-(0.12*exp(-3*Bi(i))));
  tum2=tum2+(x(i)*m(i)*(d(i)^2));
  tum3=tum3+(x(i)*m(i)*(d(i)^3)); 
end
zeta2=nu*coef*tum2;
zeta3=nu*coef*tum3;
for i=1:1:n
 for j=1:1:n   
   eps(i,j)=(epsilon(i)+epsilon(j))/2;
   be(i,j)=(sqrt(sigma(i)*sigma(j))/(0.5*(sigma(i)+sigma(j))));
   kap(i,j)=(sqrt(kappa(i)*kappa(j)))*(be(i,j)^3);
   d(i,j)=0.5*(d(i)+d(j));
   g(i,j)=(1/(1-zeta3))+(((3*d(i)*d(j))/(d(i)+d(j)))*(zeta2/(1-zeta3)^2))+(((d(i)*d(j))/(d(i)+d(j)))^2)*((zeta2^2)/(1-zeta3)^3);
 end  
end
delta11=g(1,1)*((exp(eps(1,1)/(1*T)))-1)*kap(1,1)*(d(1,1)^3);
delta12=g(1,2)*((exp(eps(1,2)/(1*T)))-1)*kap(1,2)*(d(1,2)^3);
delta22=g(2,2)*((exp(eps(2,2)/(1*T)))-1)*kap(2,2)*(d(2,2)^3);
moadele=[coe(1)+(nu*Na*((x(1)*coe(1)*coe(2)*delta11)+(x(2)*coe(1)*coe(4)*delta12)))-1
         coe(2)+(nu*Na*((x(1)*coe(1)*coe(2)*delta11)+(x(2)*coe(2)*coe(3)*delta12)))-1
         coe(3)+(nu*Na*((x(1)*coe(2)*coe(3)*delta12)+(x(2)*coe(3)*coe(4)*delta22)))-1
         coe(4)+(nu*Na*((x(1)*coe(1)*coe(4)*delta12)+(x(2)*coe(3)*coe(4)*delta22)))-1];
