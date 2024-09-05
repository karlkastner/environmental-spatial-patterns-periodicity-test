% 2023-05-16 17:34:27.559560647 +0200
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
%% demonstration, that the periodicity test is able to detect periodic components
%% in stochastic patterns
%% 
%% function example_detection_of_periodic_components_in_stochastic_patterns(pflag)
function example_detection_of_periodic_components_in_stochastic_patterns(pflag)

	if (nargin()<1)
		pflag = false;
	end
	fflag = pflag;

	% characteristic wavelength
	lambda_c = 1;
	% characteristic frequency	
	fc = 1/lambda_c;
	% regularity
	reg = 1;
	% domain size
	L = 10;
	% sampling interval
	dx = lambda_c/100;
	% number of grid points
	n = L/dx;
	% approximate fraction of spectral energy contained contributed by
	% periodic components in the compbined pattern
	qS = 0.0001;
	% weight of the stochastic pattern in the combined pattern,
	% since the spectral components scales approximately as S ~ b^2
	qb = sqrt(qS);

	% reset random number gnerator, to aways produce the identical pattern
	rng(0);
	
	fx = fourier_axis(L,n);
	x=fftshift(fourier_axis(fx));
	n = length(fx);
	
	% synthesize a periodic pattern, we rotate the pattern, to make
	% it visually distinguishable from the stochastic pattern
	a = pi/3;
	% periodic pattern
	bp = sin(4*pi*(cos(a)*x+sin(a)*x')*fc);
	% threshold
	bp = bp>0;

	% generate a stochastic striped pattern with log-normal density
	% in direction perpendicular to stripes and exponential density
	% in the direction perpendicular to the stripes
	[a,b] = lognmirroredpdf_mode2par(fc,0.5*reg/fc);
	Sx = lognmirroredpdf(abs(fx),a,b);
	[fy0,sy] = laplacepdf_mode2par(0,reg/fc);
	Sy = laplacepdf(fx,fy0,sy); % fc/4);
	S  = (cvec(Sx)*rvec(Sy));
	% uncorrelared noise
	e = randn(n);
	% transfer function of pattern
	T = sqrt(S);
	% generate the stochastic pattern by convolving the noise with
	% autocorrelation function of the process in the frequency domain
	bs = ifft2(T.*fft2(e));
	% threshold the pattern
	bs = 1-(bs>quantile(bs,0.66,'all'));
	% combined patterns, let the periodic part contribute roughly 5%
	
	% subtract mean	
	bs = bs-mean(bs,'all');
	bp = bp-mean(bp,'all');

	% normalize
	bs = bs/rms(bs,'all');
	bp = bp/rms(bp,'all');

	bc = (1-qb)*bs + qb*bp;

	% smoothing radius for the test
	nf = 4; %round(L/lambda_c);
	% select all frequency components for testing here
        fmsk=true(n);
	[~,pn,stat]= periodogram_test_periodicity_2d(bs, [L,L], nf, [], fmsk);
	printf('p-value of the stochastic pattern: %0.2g\n',pn);
	printf('Spectral energy contained in significant components %0.2g%%\n',stat.intShat_sig*100);
	[~,pn,stat] = periodogram_test_periodicity_2d(bp, [L,L],nf, [], fmsk);
	printf('p-value of the periodic pattern:   %g\n',pn);
	printf('Spectral energy contained in significant components %0.2g%%\n',stat.intShat_sig*100);
	[~,pn,stat] = periodogram_test_periodicity_2d(bc, [L,L],nf, [], fmsk);
	printf('p-value of the cobined pattern:    %g\n',pn);
	printf('Spectral energy contained in significant components %0.2g%%\n',stat.intShat_sig*100);

	% plot the stochastic pattern	
	splitfigure([2,2],[1,1],fflag);
	imagesc(bs);
	axis square;
	colormap(colormap_vegetation());
	axis square;
	axis off

	% plot the periodic pattern	
	splitfigure([2,2],[1,2],fflag);
	imagesc(bp);
	axis square;
	axis off
	colormap(colormap_vegetation());
	  
	splitfigure([2,2],[1,3],fflag);
	% for illustration, we add a larger periodic part, as it is 
	% otherwise imperceivable to the human eye
	imagesc(bs + 0.5*bp);
	axis square;
	axis off
	colormap(colormap_vegetation(255));
	
	%rms(abs(fft(0.05*bp(:))).^2,'all')./rms(abs(fft(b(:))).^2,'all')
	
	if (pflag)
		ps = 3.5;
		pdfprint(11,'img/pattern-random.pdf',ps);
		pdfprint(12,'img/pattern-periodic.pdf',ps);
		pdfprint(13,'img/pattern-combined.pdf',ps);
	end
end % function 
	
