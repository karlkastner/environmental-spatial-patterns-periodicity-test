% Thu  6 Jul 13:19:52 CEST 2023
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
% minimum working example

if (0)
	% this plantation pattern has significant frequency components
	filename = 'input/plantation-striped-spain_-2.685893_37.565588.png';
	% this pattern spans the full image has no mask
	mskname ='';
	f_hp = 0.01;
	f_lp = 0.1;
	axlim = [0 0.15,0,1];
	nf = 4;
else
	filename = 'input/pattern_2d_+07.79111_+047.92607.png';
	mskname  = 'input/pattern_2d_+07.79111_+047.92607_mask.tif';
	f_hp = 0.006;
	f_lp = 0.04;
	fc = 0.01;
	axlim = [0, 5,0,100];
	nf = 6.3;
end

% The simplest way to apply the test is to use the class Spatial_Pattern which
% takes care of all necessary processing step
sp = Spatial_Pattern();
sp.imread(filename);
sp.analyze_grid();
printf('p-value: %g\n',sp.stat.p_periodic())

% alternatively, the steps can be performed manually:

% load an image with pattern
img = imread(filename);
if (~isempty(mskname))
	bmsk = single(imread(mskname));
else
	bmsk = [];
end
% convert to grayscale,
% if a nir-band is available, use NDVI here
b = single(im2gray(img));

% the spurious low frequency lobe has to be excluded, this can be determined
% by visually inspecting the radial periodogram 
n = size(b);
% dummy image size, the true size is not required
L       = n;
mu_b    = sum(b.*bmsk,'all')/sum(bmsk,'all');
Shat    = periodogram_2d(bmsk.*(b-mu_b));
[Sr,fr,Cr] = periodogram_radial(Shat,L);

subplot(2,2,2)
plot(fr/fc,Sr.normalized)
axis(axlim);
[fx,fy,frr] = fourier_axis_2d(L,n);
% We cut the lobe at the frequency f_hp, the minimum between
% the low frequency lobe and the main lobe of the characteristic frequency
% component:
% cutting high frequency components makes the test more significant
vline(f_hp/fc,'color','r');
vline(f_lp/fc,'color','b');
yyaxis right
plot(fr/fc,Cr);
xlabel('Wavenumber k/k_c');
ylabel('Density S_r/\lambda_c');

% exclude the lobe
fmsk = frr<f_hp;
% exclude the symmetric half plane
fmsk = fmsk & (fx >= 0);
% optionally exclude high frequency components
fmsk = fmsk & frr < f_lp;

% selected a smoothing radius for estimating the density,
% we choose nf ~ sqrt(2 fc/df)
[~, pn, stat] = periodogram_test_periodicity_2d(b, [], nf,bmsk,fmsk);
% note that the p-value slightly differs from that computed by analyze_grid,
% due to the slighlty different automatic choices for nf, f_hp and f_lp
printf('p-value: %g\n', pn);

