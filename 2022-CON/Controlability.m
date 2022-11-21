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
   x2*Re/x1*(-Re*x2*mdote)/(vcc*(x3-Re))];

g=[Re/(vcc*(x3-Re))*(cpa*Ta) , Re/(vcc*(x3-Re))*(hf);  
   x2*Re/x1*( ((Re-x3)*x2+cpa*Ta))/(vcc*(x3-Re)) ,  x2*Re/x1*(((Re-x3)*x2+hf))/(vcc*(x3-Re))];
%%
%Equation 32 - C=[ g1, g2, [g1,g2], [f,g1], [f,g2], [f,[f,g1]], [f,[f,g2]]  ]

g1g2brac= [diff(g(:,2),x1),diff(g(:,2),x2)]*g(:,1)- [diff(g(:,1),x1),diff(g(:,1),x2)]*g(:,2); %[g1,g2]
fg1brac = [diff(g(:,1),x1),diff(g(:,1),x2)]*f     - [diff(f,x1),diff(f,x2)]*g(:,1);           %[f, g1]
fg2brac = [diff(g(:,2),x1),diff(g(:,2),x2)]*f     - [diff(f,x1),diff(f,x2)]*g(:,2);           %[f, g2]
ffg1brac= [diff(fg1brac,x1),diff(fg1brac,x2)]*f   - [diff(f,x1),diff(f,x2)]*fg1brac;          %[f,[f,g1]]
ffg2brac= [diff(fg2brac,x1),diff(fg2brac,x2)]*f   - [diff(f,x1),diff(f,x2)]*fg2brac;          %[f,[f,g2]]

C=[g(:,1) g(:,2) g1g2brac fg1brac fg2brac ffg1brac ffg2brac];
%%
%%%%%%%%%%%%%%%%%%%%%%%
%% Condition of the  Controllability
%%%%%%%%%%%%%%%%%%%%%%%
n=2;% number of the states of the combustion chamber 
if n==rank(C)
    fprintf('The gas turbine is  locally controllable\n');
else
    fprintf('The gas turbine isnot locally controllable\n');
end
%%%%%%%%%%%%%%%%%%%%%%%