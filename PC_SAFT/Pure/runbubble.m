global kij
global T
global x
T=313.15;
x(1)=0.1905;
x(2)=1-x(1);
kij=0.0085;
coe=[10 0.9 0.1];
jav=fsolve(@f,coe)
