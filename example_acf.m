% 2024-02-07 14:24:35.275866075 +0100
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

	if (~exist('pflag'))
		pflag = 0;
	end
	fflag = pflag;
	% characteristic wavelength
	lc=1;
	% characteristic frequency
	fc=1/lc;
	% spatial extent
	L = 40*lc;
	% number of grid points
	n=L^2;
	% generate an iostropic pattern
	[z,x,y,Lx,Ly]=generate_isotropic_pattern(1./lc,n,L);
	% periodogram
	hatS =abs(fft2(z-mean(z,'all'))).^2;
	% autocorrelation
	R = ifft2(hatS);
	% radial autocorrelation
	[Rr,r] = autocorr_radial(R,[Lx,Ly]);
	% radial periodogram
	[Sr,fr] = periodogram_radial(hatS,[Lx,Ly]);

	% plot the pattern
	splitfigure([2,3],[1,1],fflag);
	imagesc(y,x,z);
	xlabel('x/\lambda_c');
	ylabel('y/\lambda_c');
	s = 2.5;
	axis equal
	xlim(s*[-1,1])
	ylim(s*[-1,1])
	colormap(flipud(gray))

	% plot the acf
	splitfigure([2,3],[1,2],fflag);
	sa = pi*sqrt(r/lc);
	plot(r,[Rr,Rr.*sa]);
	xlim([0,2.5]);
	legend('$R_r$','$\displaystyle \pi \sqrt{\frac{r}{\lambda_c}} R_r$','interpreter','latex');
	xlabel('Lag distance');
	
	% plot the density
	 splitfigure([2,3],[1,3],fflag);
	% plot(S1.normalized)
	stem(fr/fc,Sr.normalized/L,'.','markersize',10);
	xlim([0,2.5])
	ylabel('Density $S_r/L$','interpreter','latex');
	xlabel('Wavenumber $k_r/k_c$','interpreter','latex');


	% ringed pattern
	[z,x,y,Lx,Ly]=generate_isotropic_pattern(1./lc,n,L);
	z = z-z.^3;
	splitfigure([2,3],[1,4],fflag);
	imagesc(y,x,z);
	xlabel('x/\lambda_c');
	ylabel('y/\lambda_c');
	s = 2.5;
	axis equal
	xlim(s*[-1,1])
	ylim(s*[-1,1])
	colormap(flipud(gray))
	% periodogram
	hatS = abs(fft2(z-mean(z,'all'))).^2;	
	% autocorrelation
	R = ifft2(hatS);
	R = R/R(1);
	% radial averages
	[Rr_,r] = autocorr_radial(R,[Lx,Ly]);
	[Sr_,fr] = periodogram_radial(hatS,[Lx,Ly]);

	splitfigure([2,3],[1,5],fflag);
	plot(r/lc,[Rr,Rr_]);
	ylabel('Autocorrelation R_r');
	xlabel('Lag distance x/\lambda_c')
	xlim([0,2.5]);

	splitfigure([2,3],[1,6],fflag);
	stem(fr/fc,[Sr.normalized,Sr_.normalized]/L,'.','markersize',10);
	xlim([0,2.5])
	ylabel('Density $S_r/L$','interpreter','latex');
	xlabel('Wavenumber $k_r/k_c$','interpreter','latex');

	if (pflag)
		ps = 4;
		pdfprint(11,'img/hexagonal-pattern.pdf',ps);
		pdfprint(12,'img/hexagonal-radial-autocorrelation.pdf',ps);
		pdfprint(13,'img/hexagonal-radial-density.pdf',ps);
		pdfprint(14,'img/fairy-circle-pattern.pdf',ps);
		pdfprint(15,'img/fairy-circle-radial-autocorrelation.pdf',ps);
		pdfprint(16,'img/fairy-circle-radial-density.pdf',ps);
	end

