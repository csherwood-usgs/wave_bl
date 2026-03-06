function uw = ustwv(z,zf,Kf,omega,uinf)
% USTWV - Solve Z equation for waves
% uw = ustwv(z,zf,Kf,omega,uinf)
%
% Solve iwF = d/dz(K(dF/dz))
% Form tridiagonal coefficient matrix. b,d,a = before, diagonal, after

% csherwood@usgs.gov
mz = length(z);
dz = diff(z);
dzm = diff(zf);
dzu = dz(2:end,1);
dzl = dz(1:end-1,1);
Ku = Kf(2:end);
Kl = Kf(1:end-1);
b = [0; Kl./(dzm.*dzl); 0];
a = [0; Ku./(dzm.*dzu); 0];
d = [1;...
     complex( ( -1 ./dzm).*(...
			 Ku./dzu + Kl./dzl),...
	      -omega*ones(mz-2,1) );...
     1];
c = zeros(size(d));
c(1)= 1;
cc = sy(1,mz,b,d,a,c);
uw = uinf*real(1-cc);
return
