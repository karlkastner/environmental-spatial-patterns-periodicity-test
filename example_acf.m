	if (~exist('pflag'))
		pflag = 0;
	end
	fflag = pflag;
	lc=1;
	fc=1/lc;
	L=40;
	 n=L^2;
	 [z,x,y,Lx,Ly]=generate_isotropic_pattern(1./lc,n,L);
	 imagesc(z);
	 axis equal;
	 axis tight;
	 S =abs(fft2(z-mean(z,'all'))).^2;
	 R = ifft2(S);
	 [R1,r] = autocorr_radial(R,[Lx,Ly]);
	 S1 = periodogram_radial(S,[Lx,Ly]);
	 splitfigure([2,3],[1,1],fflag);
	 imagesc(y,x,z);
	xlabel('x/\lambda_c');
	ylabel('y/\lambda_c');
	s = 2.5;
	xlim(s*[-1,1])
	ylim(s*[-1,1])
	axis equal
	axis square
	xlim(s*[-1,1])
	colormap(flipud(gray))

	 splitfigure([2,3],[1,2],fflag);
	 sa = pi*sqrt(r/lc);
	 plot(r,[R1,R1.*sa]);
	 xlim([0,2.5]);
	 legend('$R_r$','$\displaystyle \pi \sqrt{\frac{r}{\lambda_c}} R_r$','interpreter','latex');
	 xlabel('Lag distance');
	 splitfigure([2,3],[1,3],fflag);
	% plot(S1.normalized)
	stem(fr/fc,S1.normalized/L,'.','markersize',10);
%,S1_.normalized]/L);
	xlim([0,2.5])
	ylabel('Density $S_r/L$','interpreter','latex');
	xlabel('Wavenumber $k_r/k_c$','interpreter','latex');


if (1)	
	 [z,x,y,Lx,Ly]=generate_isotropic_pattern(1./lc,n,L);
	% z=z-z.^2;
	%z = z.^2;
	%z = sqrt(z)-z.^2;
	z = z-z.^3;
	 splitfigure([2,3],[1,4],fflag);
	 imagesc(y,x,z);
	xlabel('x/\lambda_c');
	ylabel('y/\lambda_c');
	s = 2.5;
	xlim(s*[-1,1])
	ylim(s*[-1,1])
	axis equal
	axis square
	xlim(s*[-1,1])
	colormap(flipud(gray))
	S = abs(fft2(z-mean(z,'all'))).^2;	
	R = ifft2(S);
	R = R/R(1);
	[R1_,r] = autocorr_radial(R,[Lx,Ly]);
	[S1_,fr] = periodogram_radial(S,[Lx,Ly]);
	splitfigure([2,3],[1,5],fflag);
	plot(r/lc,[R1,R1_]);
	ylabel('Autocorrelation R_r');
	xlabel('Lag distance x/\lambda_c')
	xlim([0,2.5]);
	splitfigure([2,3],[1,6],fflag);
	stem(fr/fc,[S1.normalized,S1_.normalized]/L,'.','markersize',10);
	xlim([0,2.5])
	ylabel('Density $S_r/L$','interpreter','latex');
	xlabel('Wavenumber $k_r/k_c$','interpreter','latex');
end
if (0)	
	L = 10*lc; 
	%ll = [lc,lc+0.5*lc/L];
	x = (0:n-1)'/n*L;
	%z = 0.5*(1+sin(2*pi*x./lc)).^3; % + 0.25*sin(2*pi*x);
	%z = [(sin(2*pi*x./lc)) -1/2*sin(4*pi*x/lc)+1/3*sin(6*pi*x/lc)];
	z = 0*(1+cos(2*pi*x)) + -(1+cos(2*pi*x)).^3;
	% 1 + 2 cos + cos^2 = 1 + 2 cos + 1/2 + 1/2 cos 2x
	%k = 1:10;
	%z = sum((-1).^k.*sin(2*pi*x./lc*k)./k,2);
	%z = sin(2*pi*x) - 0.5*(cos(4*pi*x)) + 1/3*cos(6*pi*x);
	%z = sin(2*pi*x) - 0.5/2*(cos(4*pi*x)) + 1/3/3*cos(6*pi*x);
	w = tukeywin(n);
	wz =z-mean(z);
	wz =w.*wz;
	S = abs(fft(z-mean(z))).^2;	
	R = ifft(S);
	R = R/R(1);
	 splitfigure([2,3],[1,4],fflag);
	plot(x,z)
	 splitfigure([2,3],[1,5],fflag);
	plot(x./lc,R);
%	plot(x./ll,z);
	xlim([0,2.5]);
	 splitfigure([2,3],[1,6],fflag);
	plot(S)
end
if (0)	
	L = 10*lc; 
	ll = [lc,lc+0.5*lc/L];
	x = (0:n-1)'/n*L;
	z = sin(2*pi*x./ll);
	w = tukeywin(n);
	wz =z-mean(z)
	wz =w.*wz;
	S = abs(fft(wz-mean(wz))).^2;	
	R = ifft(S);
	R = R/R(1);
	 splitfigure([2,3],[1,5],fflag);
	plot(x./ll,R);
%	plot(x./ll,z);
	xlim([0,2.5]);
	 splitfigure([2,3],[1,6],fflag);
	plot(f/fc,S)
end	

	if (pflag)
		ps = 4;
		pdfprint(11,'img/hexagonal-pattern.pdf',ps);
		pdfprint(12,'img/hexagonal-radial-autocorrelation.pdf',ps);
		pdfprint(13,'img/hexagonal-radial-density.pdf',ps);
		pdfprint(14,'img/fairy-circle-pattern.pdf',ps);
		pdfprint(15,'img/fairy-circle-radial-autocorrelation.pdf',ps);
		pdfprint(16,'img/fairy-circle-radial-density.pdf',ps);
	end

