% 2023-05-06 16:22:14.646931005 +0200
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
%% demonstration
%%
%% 1) that the critical values for the confidence levels a1 suggested by
%%    Renshaw and Ford 1984 result in the misidentification of significant frequency
%%    components, in nearly all patterns consisting of white noise,
%%    i.e. P(p>a1) ~ 1 >> a1
%%
%% 2) that the amended confidence interval with level an = a1/nt,
%%    correctly results in the identification of significant components in only
%%    P(p>an) = a1 of all patterns

% number of cells of the grid
nn = 32^2;

% numer of samples for the mc-simulation
nmc=1e4;

% significance level for testing 1 frequency component
p1 = 0.05;

% white noise "pattern"
x = randn(nn,nmc);

% periodogram
Shat = abs(fft(x)).^2;

% normalized periodogram
Shat = 2*Shat./sum(Shat);

% critical value suggested by Renshaw and Ford 1984
q_rf = 0.59/100;
p_rf = 1-chi2cdf((q_rf*nn),2);

% significance level corrected for multiple testing
pn = 1-(1-p1)^(1/(0.5*nn));
pn_ = p1/(0.5*nn);
% critical value corrected for multiple testing all nn/2 components
q_theory  = 0.5*chi2inv(1-pn,2)/(0.5*nn);

% empirical confidence level for testing all nn/2 components
q_mc = quantile(max(Shat),1-p1);

pt_theory = mean(max(Shat>q_theory));
pt_mc     = mean(max(Shat>q_mc));
pt_rf     = mean(max(Shat>q_rf));

printf('Fraction of patterns with significant components at level %g:\n',p1);
printf('Expected:                                            %g%%\n',100*p1);
printf('Estimated by Monte-Carlo simulation:\n');
printf('With critical values corrected for multiple testing: %g%%\n',100*pt_mc);
printf('With critical values suggested by Renshaw and Ford : %g%%\n',100*pt_rf);


