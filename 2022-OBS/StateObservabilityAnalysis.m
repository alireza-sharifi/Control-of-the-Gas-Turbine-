clear all
%%
%Define Variables
syms Pn x3 cpa alphap beta hf T1 T2 pamb vcc Ra Re etac etat m0 TetaccN pccN Tetaamb
syms x1 x2 x3
%%
mdote=m0*(sqrt(TetaccN)/sqrt(x2))*(x1/pccN);
Ta=   Tetaamb/etac*((x1/pamb)^(Ra/cpa)+etac-1);
Tex=  x2*etat*((pamb/x1)^(Re/x3)+1/etat-1);
%%
%Equation 10 - xdot=f(x)+g(x)u
f=[Re/(vcc*(x3-Re))*(-x3*x2*mdote);
   x2*Re/x1*(-Re*x2*mdote)/(vcc*(x3-Re));
   0];

g=[Re/(vcc*(x3-Re))*(cpa*Ta) , Re/(vcc*(x3-Re))*(hf);  
   x2*Re/x1*( ((Re-x3)*x2+cpa*Ta))/(vcc*(x3-Re)) ,  x2*Re/x1*(((Re-x3)*x2+hf))/(vcc*(x3-Re));
   0,0];
%%
%Equation 11 - z=h(x)
h=[alphap*(x1*x2);
   Tex];

Lfh1=[diff(h(1),x1) ,diff(h(1),x2),diff(h(1),x3)]*f;
Lfh2=[diff(h(2),x1) ,diff(h(2),x2),diff(h(2),x3)]*f;
Lf2h1=[diff(Lfh1,x1),diff(Lfh1,x2),diff(Lfh1,x3)]*f;
Lf2h2=[diff(Lfh2,x1),diff(Lfh2,x2),diff(Lfh2,x3)]*f;
Lf3h1=[diff(Lf2h1,x1),diff(Lf2h1,x2),diff(Lf2h1,x3)]*f;
Lf3h2=[diff(Lf2h2,x1),diff(Lf2h2,x2),diff(Lf2h2,x3)]*f;
L1gh1=[diff(h(1),x1) ,diff(h(1),x2),diff(h(1),x3)]*g;
L1gh2=[diff(h(2),x1) ,diff(h(2),x2),diff(h(2),x3)]*g;
L1fL1gh1=[[diff(L1gh1(1),x1) ,diff(L1gh1(1),x2) ,diff(L1gh1(1),x3)]*f    ,    [diff(L1gh1(2),x1) ,diff(L1gh1(2),x2) ,diff(L1gh1(2),x3)]*f];
L1fL1gh2=[[diff(L1gh2(1),x1) ,diff(L1gh2(1),x2) ,diff(L1gh2(1),x3)]*f    ,    [diff(L1gh2(2),x1) ,diff(L1gh2(2),x2) ,diff(L1gh2(2),x3)]*f];
L1gL1fh1=[diff(Lfh1,x1) ,diff(Lfh1,x2),diff(Lfh1,x3)]*g;
L1gL1fh2=[diff(Lfh2,x1) ,diff(Lfh2,x2),diff(Lfh2,x3)]*g;
%%
% Equation 16 
% W=[  0,0                       
%      L1gh                         
%     (L1gLf1 + L1fL1g)h        
%     (L1gLf2 + L2fL1g)h  ]; 
W=[  0,0                   ;
     0,0                   ;
     L1gh1,                ;
     L1gh2,                ;
     L1gL1fh1 + L1fL1gh1  ;     
     L1gL1fh2 + L1fL1gh2  ];  
%% 
% Equation 10
% I=[h(1) , h(2);
%    Lfh1 , Lfh2;
%    Lf2h1, Lf2h2];

I=[h(1) , h(2);
   Lfh1 , Lfh2;
   Lf2h1, Lf2h2];

O=[diff(I(:,1),x1), diff(I(:,1),x2), diff(I(:,1),x3);
   diff(I(:,2),x1), diff(I(:,2),x2), diff(I(:,2),x3)];
%%
%% Condition of the State Local Observability Analysis
%%%%%%%%%%%%%%%%%%%%%%%
n=3;%number of states and the gas turbine
if   (rank(O)==n && n==(rank([O  W])-rank(W)))
    fprintf('The gas turbine estimation problem is locally observable\n');
else
    fprintf('The gas turbine estimation problem isnot locally observable\n');
end
%%%%%%%%%%%%%%%%%%%%%%%

