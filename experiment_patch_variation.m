% 2023-11-30 15:21:15.015367143 +0100
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


% variation of location
% variation of size
% variation of shape (not possible) in 1d
% addition of  noise
if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;
% reset random generator
rng(0);
for jdx=1:5

n = 1e4;
m = 1e2;
% spatial resolution
dx = 1;
% characteristic wavelength
lc = m;
% characteristic frequency
fc = 1./lc;
% spatial extent
L = n*m;

b = zeros((n+2)*m,1);
M = [];
xi = 0;
 for idx=1:n-2
	switch(jdx)
	case {1}
		% fully deterministic
		b = b(:);
		k = (idx-1)*m+1;
		k_ = (0:m-1)';
		fdx = k+k_>0;
		k_ = k_(fdx);
		b(k+k_) = b(k+k_) + betapdf(k_/(m-1),5,5)/betapdf(0.5,5,5);
	case {2}
		% vary width
		m_ = round(m*(0.5+betarnd(5,5)));
		M(idx) = m_;
		k = (idx-0.5)*m+round(m_/2);
		k_ = (0:m_-1)';
		fdx = k+k_>0;
		k_ = k_(fdx);
		b(k+k_) = b(k+k_) + betapdf(k_/(m_-1),5,5)/betapdf(0.5,5,5);
	case {3}
		% vary location
		k = 1 + (idx-1)*m + round(m*(betarnd(5,5)-0.5));
		k_ = (0:m-1)';
		fdx = k+k_>0;
		k_ = k_(fdx);
		b(k+k_) = b(k+k_) + betapdf(k_/(m-1),5,5)/betapdf(0.5,5,5);
	case {4} % add noise
		b = b(:);
		k = (idx-1)*m+1;
		k_ = (0:m-1)';
		fdx = k+k_>0;
		k_ = k_(fdx);
		%b(k+k_) = b(k+k_) + betapdf(k_/(m-1),5,5)/betapdf(0.5,5,5).*(1+2*(betarnd(5,5,m,1)-0.5)); 
		b(k+k_) = b(k+k_) + betapdf(k_/(m-1),5,5)/betapdf(0.5,5,5)+(0.33*(betarnd(5,5,m,1)-0.5)/beta_std(5,5)); 
	case {5} % fully stochastic
		p = 0.0;
		q = 0.0;
		xi = xi + m*(p + (1-p)*(0.5+betarnd(2,2)));
		mi = round(m*(q + (1-q)*(0.5+betarnd(2,2))));
		k  = round(xi/dx);
		k_ = (0:mi-1)';
		x_ = k_/(mi-1);
		b(k+k_) = b(k+k_) + betapdf(x_,5,5)/betapdf(0.5,5,5);
	end
end

b = b(m+1:(n+1)*m);
% bar = mean(b);
bar = 1;
x = (0:n*m-1)';

% estiamte density with bartlett's method
if (1)
% split the pattern into sqrt-n slices
b_ = reshape(b,[],sqrt(n));
% spatial coordinate axis for a single slice
x_ = x(1:sqrt(n)*m);
% fourier axis for a single slice
fx_ = fourier_axis(x_);
% spectral resolution (equal to 1/L_)
dfx_ = fx_(2)-fx_(1);
% length of a slice
L_ = x_(end);
% periodogram of slices
hatS = abs(fft(b_-mean(b_,1))).^2;
% density estimated as average of slices
barS = mean(hatS,2);
% normalize the debsity
barS = 2*barS/(sum(barS)*dfx_);
% autocorrelation function
R = ifft(barS);
R = R/R(1);
% density maximum
[Sc,mdx] = max(barS);
% characteristic frequency
fc = abs(fx_(mdx));
% characteristic wavenumber
lc = 1./fc;
else
	 S = abs(fft(b-mean(b))).^2;
	 df = 1./L;
	 fx = fourier_axis(L,n*m);
if (jdx == 5)
	S(abs(fx)<0.025*fc) = 0;
end
	R =ifft(S);
	R = R/R(1);
if (jdx == 5)
	S=ifftshift(trifilt1(fftshift(S),200+1));
	[Sc,mdx] = max(S);
	fc = abs(fx(mdx));
	lc = 1./fc;
end
	S = 2*S/(sum(S)*df);
end

splitfigure([5,4],[1,1+4*(jdx-1)],fflag,[],100);
a=area(x(1:12*m)/lc,b(1:12*m).*[0,1]/bar);
a(2).FaceColor = [0,0.66,0];
ylim([-0.001,1.4])
ylabel('Pattern b');
xlabel('Distance x/\lambda_c');
xlim([0,8]);
axis square

splitfigure([5,4],[1,2+4*(jdx-1)],fflag,[],100);
if (5 ~= jdx)
	plot(fx_/fc,barS/L_,'linewidth',1);
	ylabel('Density S/L');
	ylim([0,1.05])
else
	yyaxis left
	%s = 1e4;
	plot(fx_/fc,barS/L,'linewidth',1);
	ylabel('Density S/L');
	yl = [0,2]*1e-4;
	ylim(yl);
	yyaxis right
	ylabel('Density S/\lambda_c');
	ylim(yl*L/lc);
	set(gca,'ycolor','k')
	% \cdot 10^{-4}');
if (0) 
	plot(fx_/fc,barS/lc,'linewidth',1);
	yyaxis right;
	ylim(yl);
	%ylim(yl*lc/L*1e-4);
	% \cdot 10^{-4}');
	set(gca,'ycolor','k')
end
	%ylim([0,1.05])
end
xlim([0 4.25]);
xlabel('Wavenumber k/k_c');
axis square

splitfigure([5,4],[1,3+4*(jdx-1)],fflag,[],100);
cla
plot(x_/lc,R,'linewidth',1);
ylim([-1.1,1.1])
xlim([0,4.25]);
hline(0)
xlabel('x/\lambda_c');
ylabel('Autocorrelation R');
axis square
hold on;
xi = linspace(0.75,5);
for idx=1:5
	 fdx = find(x_>(idx-0.5)*lc & x_<(idx+0.5)*lc);
	 [Rc(idx),mdx] = max(R(fdx));
	 xc(idx) = x_(fdx(mdx));
end
Ri = interp1(xc/lc,Rc,xi,'spline');
plot(xi,Ri,'-r','linewidth',1)

splitfigure([5,4],[1,4+4*(jdx-1)],fflag,[],100);
p = patch_size_1d(b>mean(b));
plot((1:m)/lc,lc*p(1:m)/sum(p))
axis square

end

if (pflag)
	ps = 3.5
	C = {'pattern','density','acf'};
	for idx=1:3
	pdfprint(101+idx-1,['img/experiment-1d-variation-periodic-',C{idx},'.pdf'],ps);
	pdfprint(105+idx-1,['img/experiment-1d-variation-patch-width-',C{idx},'.pdf'],ps);
	pdfprint(109+idx-1,['img/experiment-1d-variation-patch-spacing-',C{idx},'.pdf'],ps);
	pdfprint(113+idx-1,['img/experiment-1d-variation-noise-',C{idx},'.pdf'],ps);
	pdfprint(117+idx-1,['img/experiment-1d-variation-stochastic-',C{idx},'.pdf'],ps);
	end
end


