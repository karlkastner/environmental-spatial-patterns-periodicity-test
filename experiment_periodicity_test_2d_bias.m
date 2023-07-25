% Mon 15 May 11:30:20 CEST 2023
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
%
%% demonstration of the bias of the periodicity test
%
if (~exist('pflag','var'))
	pflag = false;
end
fflag = pflag;
%m=1e3;
m = 1e2;
% characteristic wavelength
lc = 100;
% characteristic frequency
fc = 1/lc;
% spatial extend of pattern, experiment is repeated for 3 length
LL = [5,10,20].*lc;
% spatial discretization
dx = 4;
% regularity of pattern, experiment is repeated for a range of regularities
reg = 2.^(-2:0.5:1);
significance_level = 0.05;
% the computation takes long, so only compute it once
if (~exist('p','var'))
	p    = zeros(m,length(reg),length(LL));
	% for each domain size
	for ldx=1:length(LL)
	L = LL(ldx);
	[fx,fy,fr] = fourier_axis_2d([L,L],[L,L]/dx);
	% spectral resolution
	df   = fx(2)-fx(1);
	% radius of disk for smoothing periodogram during density estimate
	% TODO make same choice as in analyze_grid
	% 0.5*f50 = 0.25*fc
	nf   = round(2*sqrt(L*fc));
	nf   = max(3,nf);
	ns   = [];
	% frequency range, tested for periodic components
	fmsk = (fx>=0) & (fr < 3*fc);
	% pattern mask, test whole pattern
	bmsk = [];
	%bmsk = true(L/dx);
	% for each frequency component
	for idx=1:length(reg)
	disp([ldx,idx]);
	% reset random number generator
	rng(0);
	% spectral density
	[a,b] = logn_mode2par(fc,reg(idx)/fc);

	S  = lognpdf(fr,a,b);
	% transfer function
	T = sqrt(S);
	% generate m-patterns and test them
	for jdx=1:m
		% uncorrelated (white) noise
		e = randn(size(fr));
		% generate pattern by convolving the noise with the autocorrelation function
		% in frequency space
		b = real(ifft2(T.*fft2(e)));
		% test for periodicity
		[issignificant,stat]= periodogram_test_periodicity_2d(b, [L,L], nf, bmsk, fmsk, ns);
	        p(jdx,idx,ldx) = stat.pn;
	end % for jdx
	% plot a pattern with reg ~ 1.4 and L = 10 for illustration
	if (abs(reg(idx) - sqrt(2)) < sqrt(eps) && LL(ldx) ==  10*lc)
		figure(1e3);
		clf();
		nx = L/dx;
		% this is off by one but ensures the 10 is plot
		x = (0:nx-1)/(nx-1)*LL(ldx).*fc;
		b = (b>quantile(b,0.8,'all'));
		imagesc(x,x,b);
		axis xy
		axis equal;
		%axis off
		xlabel('x/\lambda_c')
		ylabel('y/\lambda_c')
		colormap(colormap_vegetation)
		axis tight
		set(gca,'xtick',0:2:10,'ytick',0:2:10)
		if (pflag)
		ps = 3.5;
		pdfprint(1000,'img/synthesized-isotropic-reg-1.41.pdf',ps);
		end
	end
	end % for idx
	end % for ldx
	end % if ~exist p

% plot the fraction of patterns with significant frequency components
for ldx=1:length(LL)
 L = LL(ldx);
 splitfigure([2,2],[1,ldx],fflag);
 cla();
 [fx,fy,fr] = fourier_axis_2d([L,L],[L,L]/dx);
 df   = fx(2)-fx(1); % = 1/L
 nf   = round(0.5*sqrt(fc/df));
[a,b] = gamma_mode2par(fc,1/fc);
S = gampdf(fr,a,b);
% a = bandpass1d_continuous_pdf_max2par(fc,1/fc);
% S = bandpass1d_continuous_pdf(fr,fc,a,1);
 [Sbar,nf2] = circfilt2(S,nf);
 [Sr,fr] = periodogram_radial(S,[L,L]);
 Sbarr = periodogram_radial(Sbar,[L,L]);
 plot(fr/fc,[Sr.normalized,Sbarr.normalized]*fc,'linewidth',1);
 xlim([0,2.5]);
 xlabel('Wavenumber k / k_c');
 ylabel('Radial density S_r / \lambda_c');
 legend('$S_r$','$\bar S_r$','interpreter','latex');
 axis square
end

 splitfigure([2,2],[1,4],fflag);
 cla()
 P  = squeeze(mean(p<=significance_level));
 sP = sqrt(P.*(1-P));
 plot(reg,P,'*');
 set(gca,'xscale','log');
 d=0.25;
 xlim([min(reg)-d,max(reg)+d]);
 xlabel('Regularity S_c/\lambda_c');
 hline(significance_level,'color','k','linestyle','--','linewidth',1);
 ylim([0,0.25]);
 set(gca,'xtick',roundn(reg,2));
 ylabel('P$(p\le0.05)$','interpreter','latex');
 ylim([0,0.5]);
 xlim([0.25/2^(1/4),2*2^(1/4)]);
 hline(0.05);
 lh = legend(num2str(cvec(LL*fc)),'location','northwest');
 title(lh,'$L/\lambda_c$','interpreter','latex');
 axis square
if (pflag)
	ps = 3.5;
	pdfprint(11,'img/estimation-bias-radial-density.pdf',ps);
	pdfprint(14,'img/estimation-bias-periodicity-test.pdf',ps);
end

