clear all
%close all

% INPUT BEGIN %%%%%%%%%%%%%%%
L=1;
a=10;
kappa=1e-1;
Nel=20;
phi1=0; %BC Left
phi2=1; %BC Right
kappa_tilde_switch = 1; % -> use "1" for Galerkin+SUPG and "0" for Galerkin
% INPUT END %%%%%%%%%%%%%%%%%


PeL=a*L/kappa;
h=L/Nel;
x=L/Nel*[0:1:Nel]'; 

% output
alpha = a*h/(2*kappa);
if (alpha==0)
    xi=0;
else
    coth_alpha = (exp(alpha)+exp(-alpha))/(exp(alpha)-exp(-alpha));
    xi = coth_alpha-1/alpha;
end

DC=1/2*[- 1 1;
          -1 1];

DD=1/h*[ 1 -1;
        -1  1];

M1=zeros(Nel+1,Nel+1);
for i=1:Nel
    M1(i:i+1,i:i+1) = M1(i:i+1,i:i+1) ...
                     +a*DC+(kappa+kappa_tilde_switch*a*h/2*xi)*DD;
end

M1= M1(2:end-1,:);
f11=M1(:,1);
f12=M1(:,end);
M1= M1(:,2:end-1);
phiU = [phi1; -inv(M1)*(f11*phi1+f12*phi2); phi2];


bcf=[phi1; phi2];
bcK=[1 1; 1 exp(PeL)];
ccc = inv(bcK)*bcf;
c1=ccc(1,1);
c2=ccc(2,1);
xe=[0:0.0005:L];

figure(2), clf
gg=plot(x,phiU,'.-b');
set(gg,'LineWidth',1) 
hold on
gg=plot(xe,c1+c2*exp(PeL*xe/L),'k');
set(gg,'LineWidth',1)
axis([0 1 -.2 1])

