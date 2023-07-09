% Mon 15 May 12:56:12 CEST 2023
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
%% sketch of frequency range the periodicity test is restricted to
%
function example_frequency_range_of_test(pflag)
	if (nargin()<1)
		pflag = false;
	end

	% characteristic wavelength
	lambda_c = 1;
	% characteristic frequency
	fc=1/lambda_c;
	% domain size
	L = 10*lambda_c;
	% spatial discretization
	dx=lambda_c/100;
	% regularity
	reg = 1;
	fx = fftshift(fourier_axis(L,L/dx));
	fr=fx(fx>=0);
	% isotropic log-normal-density
	% combined with a spurious low frequency lobe in shape of an exponential 
	[a,b] = logn_mode2par(fc,reg/fc);
	q=0.5;
	% spectral density
	S = q*lognpdf(fr,a,b) + (1-q)*exppdf(fr,fc/2);
	% cumulative spectral distribution
	C = cumsum(fr.*S)*(fx(2)-fx(1));
	% normalize
	C=C/C(end);

	% lower limit
	mdx = find( S(2:end-1) < S(1:end-2) & S(2:end-1) < S(3:end),1,'first')+1;
	% upper limit
	udx = find(C>=0.8,1,'first');
	
	figure(1);
	clf();
	plot(fr/fc,S*fc,'linewidth',1);
	vline(fr(udx)/fc,'color','k','linestyle','--');
	vline(fr(mdx)/fc,'color','k','linestyle','--');
	xlim([0,2.5]);
	xlabel('Wavenumber k_r/k_c');
	ylabel('Radial density S_r \cdot f_c');
	yyaxis right
	plot(fr/fc,C,'r','linewidth',1);
	ylabel('Cumulative density C_r');
	set(gca,'ycolor','r') 
	axis square
	
	if (pflag)
	 	ps = 2.5;
		pdfprint(1,'img/periodicity-test-frequency-range.pdf',ps);
	end
	
end % function

