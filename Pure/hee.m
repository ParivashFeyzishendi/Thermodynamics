function barham=hee(coe)
global kij
global T
global x
kij=0.0085;
Texp=[303.15
303.15
303.15
303.15
303.15
303.15
303.15

];
Pexp=[1.123
2.488
4.147
6.362
9.4
12.389
15.893

];
xexp=[0.0415
0.0896
0.136
0.2019
0.2753
0.3384
0.4012

];
yexp=[0.9415
0.996
0.936
0.9019
0.9753
0.9384
0.9412
];
NP=length(Texp);
for je=1:NP
    T=Texp(je,1);
    x(1)=xexp(je,1);
    x(2)=1-x(1);
    P0=Pexp(je,1)*101325;
   b=yexp(je,1);
    coe=[P0 b 1-b];
    jav=fsolve(@f,coe);
    Pcal(je,1)=jav(1);
    ycal(je,1)=jav(2);
    er(je,1)=(abs(Pexp(je,1)-Pcal(je,1))/Pexp(je,1));
    %+(abs(yexp(je,1)-ycal(je,1))/yexp(je,1));
end
sum=0;
for i=1:NP
    sum=sum+er(i,1);
end
barham=sum/NP
Pcal
ycal
end
