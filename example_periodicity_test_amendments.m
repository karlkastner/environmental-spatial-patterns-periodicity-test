% Mon 10 Apr 15:54:31 CEST 2023
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
%% plot of an example density and periodogram
%
function example_periodicity_test_amendments(pflag)
	if (nargin()<1)
		pflag = false;
	end
	fflag=pflag;
	% print scale
	ps = 3.5;

	% repeat 100 times and check 
	nn = 100;

	% number of bins along the 1d-transect
	nx = 200;
	% flat density of white noise
	Sw = 2/nx;
	% 95%-confidence interval for testing a single, a-priori known frequency components
	p1 = 0.95;
	q1 = 0.5*Sw*chi2inv(p1,2);
	% median
	pm = 0.5;
	qm = 0.5*Sw*chi2inv(pm,2);
	% 95%-confidence interval for testing all frequency components
	%a = 1-p1;
	%pn = (1-p1)^(2/n);
	%an = 1 - (1-a)^(2/nx);
	%pn = 1-an;
	nt = nx/2;
	%an = a/(nt);
	pn = 1 - (1-p1)/nt;
	%pn = 1 - p1/(n);
	qn = 0.5*Sw*chi2inv(pn,2);
	nf = 5;

	% reset random number generator for reproducing everytime the same result
	rng(0);

	fx = fourier_axis(1,nx);
	df = fx(2)-fx(1);
	fdx  = (fx>=0);
	PP = zeros(nn,3);
	for idx=1:nn

	% pattern consisting of white noise
	b = randn(nx,1);
	
		
	% periodogram
	Shat = abs(fft(b)).^2;
	% normalize
	Shat = 2*Shat/(sum(Shat)*df);
	
	% plot first pattern
	if (1==idx)
		splitfigure([2,2],[1,1],fflag);
		cla();
		h=plot(fx(fdx),Shat(fdx)/Sw,'.');
		h.HandleVisibility='off';
		ylabel('Periodogram $\hat S/\bar S$','interpreter','latex');
		xlabel('Wavenumber k');
		col = [0.5,0.5,0.5];
		%hline(1,'color',col,'linewidth',1,'linestyle','-')
		hline(q1/Sw,'color',col,'linewidth',1,'linestyle','-')
		hline(qn/Sw,'color',col,'linewidth',1,'linestyle','-')
		%legend('Median','95 single','95 all');
		o = 0.2;
		x0 = 101;
		%text(x0,    1+o,'50%-1');
		text(x0,q1/Sw+o,'95%-1');
		text(x0,qn/Sw+o,'95%-n');
		xlim([0,100])
		ylim([0,12])
	end

	% fraction frequency components, exceeding or equal to the significance level for testing a single components
	PP(idx,1) = mean(Shat(fdx)>q1);
	% number of frequency components, exceeding or equal to the significance level for testing the all components
	PP(idx,2) = sum(Shat(fdx)>=qn);
	% number of frequency components, exceeding the median
	PP(idx,3) = mean(Shat(fdx)>qm);
	
	end % for idx
	disp([median(PP(:,1)),mean(PP(:,2)>0),median(PP(:,3))]);

	% density of a stochastic process with log-normal density
	[a,b] = lognpdf_mode2par(35,1/35);
	Ss    = lognpdf(abs(fx),a,b);
	% density of a deterministic process with 1 frequency component
	Sp    = (abs(fx) == 20) + eps*Sw;
	% normalization
	Sp = 2*Sp/(sum(Sp)*df);
	Ss = 2*Ss/(sum(Ss)*df);
	% combined density
	p  = 0.025;
	S  = (1-p)*Ss+p*Sp;
	S  = 2*S/(sum(S)*df);
	Sbar = meanfilt1(Ss,nf);

	splitfigure([2,2],[1,2],fflag);
	[hn,h] = hist(PP(:,1)*100,0.5:1:10.5);
	bar(h,hn/sum(hn));
	xlabel('Number of frequency components exeeding \alpha_1');
	
	splitfigure([2,2],[1,3],fflag);
	[hn,h]=hist(PP(:,2));
	bar(h,hn/sum(hn));
	xlabel('Number of frequency components exeeding \alpha_n');
	hline(0.05);

	splitfigure([2,2],[1,4],fflag);
	p = (1:nx/2)'/(nx/2+1);
	q = 0.5*chi2inv(p,2);
	plot([q,sort(Shat(fdx))/Sw],p);
	title('Cumulative distribution');
	legend('Theoretical \chi^2','simulated','location','southeast');
	xlabel('\hat S');
	ylabel('p');


%clf
%plot([S,Sp,Ss])
%[sum(S),sum(Sp),sum(Ss)]
%pause
	
	splitfigure([2,2],[2,1],fflag);
	plot(fx(fdx),[S(fdx), Sp(fdx),Ss(fdx)]/Sw,'linewidth',1);
	hline(1,'linewidth',1,'linestyle','--','color','k');
	set(gca,'colororder',[0,0,0; 0.8,0,0; 0,0,0.8])
	xlabel('Wavenumber k');
	ylabel('Density $S / S_{\textrm{w}}$','interpreter','latex');
	legend('Combined','Periodic','Stochastic','White')
	set(gca,'yscale','log')
	ylim([0.001,1]/Sw)
	xlim([0,nx/2])
	axis square
	xlim([0,100]);
	
	Shat = Shat.*S;
	Shat = Shat/(sum(Shat(fdx))*df);

	splitfigure([2,2],[2,2],fflag)
	cla();
	plot(fx(fdx),Shat(fdx)/Sw,'.')
	hold on
	%plot(fx(fdx),Ss(fdx)/S_,'b-','linewidth',1)
	%plot(fx(fdx),Sp(fdx)/S_,'r--','linewidth',1)
	if (1)
	%hline(1,'color',col,'linewidth',1,'linestyle','-')
	hline(q1/Sw,'color',col,'linewidth',1,'linestyle','-')
	hline(qn/Sw,'color',col,'linewidth',1,'linestyle','-')
	o=0.6;
	%text(79,    1+o,'50%-1');
	text(78,q1/Sw+o,'95%-1');
	text(78,qn/Sw+o,'95%-n');
	
	end
	%ylabel('\hat S / S_w');
	ylabel('Periodogram $\hat S / S_{\textrm{w}}$','interpreter','latex');
	xlabel('Wavenumber k');
	axis square
	xlim([0,100]);
	ylim([0,15]);
	%legend('$\hat S$ periodogram','$S$ density','interpreter','latex');
	%,'Median','95%-single','95%-all','interpreter','latex');
	
	splitfigure([2,2],[2,3],fflag)
	cla
	plot(fx(fdx),Shat(fdx)./Sbar(fdx),'.')
	hold on
	%plot(fx(fdx),Shat(fdx)./Ss(fdx),'o')
	if (1)
	%hline(1,'color',col,'linewidth',1,'linestyle','-')
	hline(q1/Sw,'color',col,'linewidth',1,'linestyle','-')
	hline(qn/Sw,'color',col,'linewidth',1,'linestyle','-')
	
	o=0.45+0.8;
	%text(79,    1+o,'50%-1');
	text(78,q1/Sw+o,'95%-1');
	text(78,qn/Sw+o,'95%-n');
	end
	xlabel('Wavenumber k');
	ylabel('$\hat S / \bar S$','interpreter','latex');
	axis square
	xlim([0,100]);
	pdx = Sw>0&fdx>0;
	%ylim([0,1.1*max(Shat(pdx)./Sbar(pdx))])
	ylim([0,15]);
	
	splitfigure([2,2],[2,4],fflag)
	cla();
	e = randn(nx,1);
	plot(real(ifft(sqrt([S,Sp,Ss]).*fft(e))))
	xlim([0,nx/6]);
	set(gca,'colororder',[0,0,0; 0.8,0,0; 0,0,0.8])

	if (pflag)
		pdfprint(11,'img/testing-levels.pdf',ps);
		pdfprint(21,'img/testing-density.pdf',ps);
		pdfprint(22,'img/testing-periodogram-raw.pdf',ps);
		pdfprint(23,'img/testing-periodogram-rescaled.pdf',ps);
	end
end % function

