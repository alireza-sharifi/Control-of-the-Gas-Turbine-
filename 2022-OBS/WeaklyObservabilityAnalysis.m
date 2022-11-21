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

g=[Re/(vcc*(x3-Re))*(cpa*Ta) , Re/(vcc*(x3-Re))*(hf),0;  
   x2*Re/x1*( ((Re-x3)*x2+cpa*Ta))/(vcc*(x3-Re)) ,  x2*Re/x1*(((Re-x3)*x2+hf))/(vcc*(x3-Re)),0;
   0,0,0];
%%
%Equation 11 - z=h(x)
h=[alphap*(x1*x2);
   Tex];
%%
% Equation 16 
% I=[  h1     h1   
%      Lfh1   Lfh2           
%      Lf2h1  Lf2h2 ]
%
% O=dI/dx 
Lfh1=[diff(h(1),x1) ,diff(h(1),x2),diff(h(1),x3)]*f;
Lfh2=[diff(h(2),x1) ,diff(h(2),x2),diff(h(2),x3)]*f;
Lf2h1=[diff(Lfh1,x1),diff(Lfh1,x2),diff(Lfh1,x3)]*f;
Lf2h2=[diff(Lfh2,x1),diff(Lfh2,x2),diff(Lfh2,x3)]*f;

I=[h(1) , h(2);
   Lfh1 , Lfh2;
   Lf2h1, Lf2h2];

% O=dI/dx
O=[diff(I(:,1),x1), diff(I(:,1),x2), diff(I(:,1),x3);
    diff(I(:,2),x1), diff(I(:,2),x2), diff(I(:,2),x3)];
%%%%%%%%%%%%%%%%%%%%%%%
%% Condition of the State Weakly Observability Analysis
%%%%%%%%%%%%%%%%%%%%%%%
n=3;% number of states and the gas turbine
if n==rank(O)
    fprintf('The states and the critical parameter of the gas turbine are weakly observable.\n');
else
    fprintf('The states and the critical parameter of the gas turbine arenot weakly observable.\n');
end
%%%%%%%%%%%%%%%%%%%%%%% 
