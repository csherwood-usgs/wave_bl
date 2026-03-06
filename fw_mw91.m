function fwc = fw_mw91( Abkn )
% fw_mw91 - Wave friction factor from Madsen and Wikramanayake (1991)

fwc = exp( 5.2*Abkn.^(-0.19)-6.1)-0.24*Abkn.^(-1.2); % Eqn. 255


if(Abkn >1000) % secant method to solve Eqn. 254 (p. 452, Numerical Recipes)
xacc = 1e-4;
MAXIT = 20;
% M = -0.1; % recommended value for MW91
M = 0.17; % matches Grant+Madsen 1979
logAbknM = log10(Abkn)-M;
fwmax = 1e-1;
fwmin = 1e-3;
x1 = 1 ./(4*sqrt(fwmax));
x2 = 1 ./(4*sqrt(fwmin));
flast = x1+log10(x1)-logAbknM;
f     = x2+log10(x2)-logAbknM;

if(abs(flast)<abs(f))
   rts = x1;
   xlast = x2;
   % swap f and flast;
   tmp = flast; 
   flast = f;
   f = tmp;
else
   xlast = x1;
   rts = x2;
end

dx = 999;
while (i<MAXIT)&(abs(dx)>xacc)
   dx = (xlast-rts)*f/(f-flast);
   xlast = rts;
   flast = f;
   rts = rts+dx;
   f = rts+log10(rts)-logAbknM;
end
fwc = (1/(4*rts)).^2;
end
