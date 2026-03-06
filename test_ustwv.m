% test_ustwv - Script to test wbl solution against analytical result
vk = 0.41;
zo = 1e-5;
zg = 2;
ustrm = .022/vk; % to match Fig. 3
w = 1.59;
uinf = 1;

% vertical grid
mz =50; % adding lots more points does not help the u*w estimates below
z = logspace( log10(zo),log10(zg), mz )';
zf = logspace( 0.5*(log10(z(1))+log10(z(2))), 0.5*(log10(z(mz))+log10(z(mz-1))), mz-1)';

% Linear eddy viscosity
Kf = vk*ustrm*zf;
% Numerical solution
uw = ustwv(z,zf,Kf,w,uinf);

% Smith (1977) solution, p. 547-548
Ko = vk*ustrm;
x = 2*sqrt(w*z/Ko);
x0 = 2*sqrt(w*zo/Ko);
% See Abramowitz & Stegun, 9.9.2, p. 379
bex = besselk(0,x.*exp(pi*1i/4));
kerx = real(bex);
keix = imag(bex);
bex0 = besselk(0,x0.*exp(pi*1i/4));
kerx0 = real(bex0);
keix0 = imag(bex0);
wt = 0;
% Eqn. 13
u2 =uinf*( cos(wt)...
   -( (kerx*kerx0+keix*keix0)/(kerx0.^2+keix0.^2) ).*cos(wt)...
   -( (kerx*keix0-keix*kerx0)/(kerx0.^2+keix0.^2) ).*sin(wt));
%% Madsen + Wikramanayake recap of GM model
del = Ko/w; % Eqn. 62
es = 2*sqrt(z./del);
es0 = 2*sqrt(zo/del);
bes = besselk(0,es.*exp(pi*1i/4));
kers = real(bes);
keis = imag(bes);
bes0 = besselk(0,es0.*exp(pi*1i/4));
kers0 = real(bes0);
keis0 = imag(bes0);
uwm = uinf*real( (1- (kers + 1i*keis)./(kers0+1i*keis0))*exp(1i*wt) );
%%
figure(1);clf
h1=semilogy(uw,z,'-k');
hold on
h2=semilogy(u2,z,'+k');
h3=semilogy(uwm,z,'ob');
legend([h1;h2;h3],'Numerical','Smith, 1977','M+W, 1991')

%%
% Ok, the velocity profiles match well, but the resulting estimates of
% bottom stress overestimate u*w by a little, even with the analytical
% solution.
ustrm
disp('All of the following values should match ustrm.')

ustrw=sqrt( Kf(1)*(abs(uw(2)-uw(1))/(z(2)-z(1))) )
ustrw2=sqrt( Kf(1)*(abs(u2(2)-u2(1))/(z(2)-z(1))) )

% Smith calculates u*w with an analytical du/dz (Eqns. 16b and 17)
k1=exp(-pi*1i/2)*besselk(1,x0*exp(pi*1i/4));
ker1 = real(k1);
kei1 = imag(k1);
F1 = ((ker1+kei1)*kerx0+(kei1-ker1)*keix0)/...
   (sqrt(2)*(kerx0.^2+keix0.^2));
F2 = (-(ker1+kei1)*keix0+(kei1-ker1)*kerx0)/...
   (sqrt(2)*(kerx0.^2+keix0.^2));
ustrm2 = sqrt(0.5*uinf*Ko*x0*sqrt(F1.^2+F2.^2))

% M+W 1991 converts equivalent of Smith's Eqn. 16 to a friction factor
Ab = uinf/w;
kb = 30*zo;
fw = fw_mw91(Ab/kb);
% This particular fw may not be an exact match, and it does not give a good
% estimate of u*w
ustrm3 = sqrt(0.5*fw*uinf.^2)


