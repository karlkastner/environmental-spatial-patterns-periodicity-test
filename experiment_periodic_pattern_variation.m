% Tue 16 Jan 10:34:11 CET 2024
% Karl Kastner, Berlin

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;
p = 0.375;
cmap = flipud((1-p)*gray()+p*colormap_vegetation(256));
% wavelength
lambda = 1;
% domain size, determines spectral resolution
L  = lambda*20*[1,1];
df = 1./L;
% grid cell size, determines spatial resolution
dx = lambda/20*[1,1];
% number of grid cells
nb = L./dx;
% spatial coordinates
x=(0:nb(1)-1)'*dx(1);
% radial coordinate
r = hypot(x-L(1)/2,x'-L(2)/2);
% frequency space coordianges f=k/2pi
fx = fourier_axis(L(1),nb(1));
ns = 100;

% plot limit
lim = 2.5;
p_thresh = [0.8,0.55];

% repeat for isotropic and anisotropic pattern
for idx=1:2
aniso = idx-1

% generate first an unperturbed and then an unperturbed image
for kdx=1:2
% reset random number generator for reprodiucibility
rng(0);

if (1 == kdx)
sn = 0;
sd = 0;
sw = 0;
sa = 0;
else
% noise level of high frequency components
sn = 0.05;
% magnitude of stripe-spacing variation
sd = 4.5;
% magnitude of random width-variation 
sw = 3.5;
sa = 0;
end
% correlation length of stripe-spacing variation
ln = 0.25*lambda;

% smoothing kernel impulse response
bs = exp(-r/ln);
bs = bs./sum(bs(:)); %*dx(1)*dx(2));
% transfer function
T  = fft2(bs);

% perturbation of x and y coordinates (phase)
% this varies the distance between patches or respectively stripes
ex = sd*randn(nb);
ex = ifft2(T.*fft2(ex)) + sn*randn(nb);
ey = sd*randn(nb);
ey = ifft2(T.*fft2(ey)) + sn*randn(nb);
% perturbartion of the stripe or respectively patch width
if (sw > 0)
[ap,bp]=logn_moment2par(2,sw);
p=lognrnd(ap,bp,nb);
p = ifft2(T.*fft2(p));
else
	p = 2*ones(nb);
end

if (aniso)
	% generate an anisotropic, i.e. striped pattern
	b = ((1+sin(2*pi/lambda*(x'+ex)))/2);
else
	% generate an isotropic, i.e. spotted pattern
	b = 0;
	xx = repmat(x,1,nb(2));
	yy = xx';
	for jdx=1:3
		R = rot2((jdx-1)*pi/3);
		xy = R*[xx(:)+ex(:),yy(:)+ey(:)]';
		x_ = reshape(xy(1,:),nb);
		b = b+cos(2*pi/lambda*x_);
	end
	b = (3+2*b)/9;
end
	% perturb band width
	%b_ = b;
	fdx = (b>1);
	b(fdx) = b(fdx).^p(fdx);
	b(~fdx) = 1 - (1-b(~fdx)).^(1./p(~fdx));
if (sa>0)
	% perturb magnitude
	[ap,bp] = logn_moment2par(1,2);
	a = lognrnd(ap,bp,nb);
	a = real(ifft2(T.*fft2(a)));
	b = b.*a;
end


	fc = 1./lambda;
	f_50 = fc;
	dfr = df(1);
	nf_test = round(0.25*f_50/dfr);
	[fx,fy,frr] = fourier_axis_2d(L,nb);
	fmsk = (frr<4*fc);
	bmsk = [];
	% TODO use Spatial_Pattern/analyse_grid here
        [isp, pn, stati] = periodogram_test_periodicity_2d(...
					b, [L,L],nf_test, bmsk, fmsk, ns);
	disp([idx,kdx,pn,stati.p1])
	stati.p1

	% threshold the pattern, this is optional
	if (1)
	q = quantile(b(:),p_thresh(idx));
	b_ = b>q;
	end
	%q = quantile(b,[0.1,1],'all');

	% periodogram, we skip the normalization here, as we do not give values for the z-scale
	Shat = abs(fft2(b-mean(b(:)))).^2;
	% normalize
	Shat = 2*Shat/sum(Shat(:)*df(1)*df(2));
	% slightly smooth it, for enlarging the tiny dots of the dominant
	% frequency components, which may vanish otherwise,
	% we compensate for the drop in magnitude which is 1/4 for nf = 3,
	% 1/9 for nf = 5
	S=9*ifftshift(trifilt2(fftshift(Shat),5));

	% corellogram (autocorrelation)
	R = ifft2(Shat);
	R = real(R);
	R = R/R(1);

	% plot the pattern
	splitfigure([2,3],[kdx,3*(idx-1)+1],fflag);
	cla();
	imagesc(x,x,b_);
	%surface(x,x,b,'edgecolor','none');
	xlabel('Distance x/\lambda_c');
	ylabel('Distance y/\lambda_c');
	axis square
	axis tight
	axis(10*[0,1,0,1]);
	axis xy
	colormap(cmap);
	%caxis(q);
	
	% plot the periodogram	
	splitfigure([2,3],[kdx,3*(idx-1)+2],fflag);
	cla();
	imagesc(fftshift(fx),fftshift(fx),fftshift(S)/lambda^2);
	axis square;
	axis tight;
	axis(lim*[-1,1,-1,1]);
	axis xy;
	colormap(cmap);
	xlabel('Wavenumber k_x/k_c');
	ylabel('Wavenumber k_y/k_c');
	c = colorbar();
	%xlim(c,[0,1]);
	set(c,'location','east')
	title(c,'$\hat S/\lambda_c^2$','interpreter','latex');

	% plot the correlogram
	splitfigure([2,3],[kdx,3*(idx-1)+3],fflag);
	cla();
	imagesc(x-L(1)/2,x-L(2)/2,fftshift(R));
	axis square;
	axis tight;
	axis(lim*[-1,1,-1,1]);
	axis xy
	colormap(cmap);
	xlabel('Distance x/\lambda_c');
	ylabel('Distance y/\lambda_c');

	c = colorbar();
	set(c,'location','east');
	title(c,'$\hat R$','interpreter','latex');
	%'xtick',[-0.5,1],'xticklabel',{'0','1'},'location','east')
	
	% k = patch_size_2d(b);
	% semilogy(trifilt1(k,111))
if (1 == kdx)
	R0 = R;
	b0 = b;
else
	disp([aniso, rms(flat(b0-b))./std(b0(:)), corr(flat(b),flat(b0)).^2])
end
end % for kdx

end % for idx  

	
if (pflag)
	ps = 3.5;
	pdfprint(21,'img/experiment-periodic-isotropic-varying-pattern.pdf',ps);
	pdfprint(22,'img/experiment-periodic-isotropic-varying-periodogram.pdf',ps);
	pdfprint(23,'img/experiment-periodic-isotropic-varying-acf.pdf',ps);
	pdfprint(24,'img/experiment-periodic-anisotropic-varying-pattern.pdf',ps);
	pdfprint(25,'img/experiment-periodic-anisotropic-varying-periodogram.pdf',ps);
	pdfprint(26,'img/experiment-periodic-anisotropic-varying-acf.pdf',ps);
end

