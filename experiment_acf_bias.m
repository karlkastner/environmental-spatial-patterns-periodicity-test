% 2024-05-28 17:37:08.844071554 +0200
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% demonstrate the bias of the acf for different methods of estimation
% when the domain length L is not an integer multiple of the characteristic
% wavenumber lambda_c

if (~exist('pflag','var'))
	pflag = 0;
end

% domain length
L   = 10;
% characteristic frequncy fc=1/lambda_c
% this is perturbed, so that fc*L = L/lc is not integer
fc0 = 1;
% perturbation of the frequency
% worst cases are +/- 1.2
m   = 100;
dfc = innerspace(-1/2,1/2,m)'/L;
% perturb frequency
fc  = fc0 + dfc;

% = [1,(L+0.5)/L];

% number of points along transect
n   = 1024*8;
% coordinates of points
x   = (0:n-1)'/n*L;
% 
fdx = x > 0.5 & x< 1.5;

a_max = [];
a = [];
% for each frequency
for idx=1:length(dfc)
	% periodic function
	y = sin(2*pi*x*fc(idx));
	% autocorrelation, estinmated in three ways
	% truncated
	a(:,1) = autocorr_unbiased(y);
	% padding zeros (matlab default)
	% this yields the same result as the matlab function in the
	% econometric toolbox
	a(:,2) = autocorr_zero(y);
	% periodically extended
	a(:,3) = autocorr_periodic(y);
	a = a./a(1,:);
if (0)
	% theoretical value
	a(:,4) = a(:,1)./(1-(1:n)'/n);
end
	% first maximum of the autocorrelation function
	a_max(idx,:) = max(a(fdx,:));
%	plot(x,a)
end

figure(1)
clf
%plot(dfc,(a_max-1)./dfc)
%plot(dfc*L/fc0,(a_max-1).*cvec(fc0+dfc)*L)
plot((L*fc - L*fc0),(a_max-1).*cvec(fc0+dfc)*L,'linewidth',1)
xlabel({'Missmatch of domain length','$L/\lambda_c - \mathrm{round}(L/\lambda_{c})$'},'interpreter','latex')
ylabel('Normalized Bias : $(R_c-1) L/\lambda_c$','interpreter','latex')
legend('Truncated (Unbiased)','Extended with mean','Extended periodically','location','south')
set(gca,'colororder',colormap_krb())
if (pflag)
	ps = 3;
	pdfprint(1,'img/acf-bias.pdf',ps);
end

