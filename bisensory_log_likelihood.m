%BISENSORY_LOG_LIKELIHOOD Calculate (log) likelihood of bimodal dataset 
% DATA under parameter set PARAMS.
%
% The values of PARAMS are:
% 1. Noise base standard deviation, visual (low noise) - log scale
% 2. Noise base standard deviation, visual (medium noise) - log scale
% 3. Noise base standard deviation, visual (high noise) - log scale
% 4. Noise base standard deviation, vestibular - log scale
% 5. Noise Weber's fraction, visual (low noise)
% 6. Noise Weber's fraction, visual (medium noise)
% 7. Noise Weber's fraction, visual (high noise)
% 8. Noise Weber's fraction, vestibular
% 9. Lapse rate (probability of random response)
% 10. Gaussian prior mean
% 11. Gaussian prior standard deviation (log scale)
% 12. Probability of common cause (`p_common`)
%
% DATA is data matrix. Each row is a trial. 
% For a given row, the columns contain data for:
% 1. Subject id (unused),
% 2. Number of trial (unused),
% 3. Task type (1 -, 2 Vestibular localisation, 3 Unity judgment)
% 4. Vestibular stimulus position (deg),
% 5. Visual stimulus position (deg),
% 6. Response (1 or -1 for right/left; 1 or 2 for yes/no unity),
% 7. Visual noise level (1 low, 2 med, 3 high).
%
% ...
function varargout = bisensory_log_likelihood(params,data,grid_size,sum_flag)

if nargin < 3 || isempty(grid_size) || isnan(grid_size); grid_size = 401; end
if nargin < 4 || isempty(sum_flag); sum_flag = true; end

bincenters = [-45,-40,-35,-30:2.5:-2.5,-1.25:0.625:1.25,2.5:2.5:30,35,40,45]';
[bimodal_counts, bincenters_bim] = data2counts(data,bincenters);
Nparams = 8; % Number of model parameters in a single noise condition

total_loglike = 0;
for iNoise = 1:3
    
    theta = zeros(1,Nparams);
    theta(1:2) = params([iNoise,iNoise + 4]);
    theta(3:4) = params([4,8]);    
    theta(5:end) = params(9:end);
    llike{iNoise} = bisensory_log_likelihood_noise( ...
        theta,bimodal_counts{iNoise},bincenters_bim,grid_size,sum_flag);
        
    if sum_flag
        total_loglike = total_loglike + llike{iNoise};
    else
        error('Individual likelihoods not supported yet.');
    end
end

varargout{1} = total_loglike;

end

%--------------------------------------------------------------------------
function varargout = bisensory_log_likelihood_noise(theta,counts,bincenters,grid_size,sum_flag)

DEBUG = 0;  % Plot some debug graphs

% Fixed parameters
MAXSD = 5;              % When integrating a Gaussian, go up to this SDs away
FIXEDLAPSEPDF = 1e-4;   % Cutoff to the penalty for a single outlier
SSCALE = 8;             % Precision of stimulus grid
MAXRNG = 90;
MAXRNG_XMEAS = 180;
alpha_rescaling_vis = 1;    % No sensory rescaling
alpha_rescaling_vest = 1;
beta_softmax = 1e4;     % Solves some numerical instabilities

% Sensory noise shape
sigma_fun = @(s, sigmazero, w) sigmazero .* (1 + (90/pi)*abs(sin(s*pi/90)).*w);

% Model parameters
lambda = theta(end-3);              % Lapse rate
prior_mu = theta(end-2);            % Gaussian prior mean
prior_sigma = exp(theta(end-1));    % Gaussian prior standard deviation
priorc1 = theta(end);               % Probability of common cause (p_common)

% Take model parameters
sigmazero_vis = exp(theta(1));
w_vis = theta(2);

sigmazero_vest = exp(theta(3));
w_vest = theta(4);

% Trials to be computed
do_estimation = ~isempty(counts{2});
do_unity = ~isempty(counts{3});

% Bin centers is a column vector
bincenters_vis = bincenters{1};
bincenters_vest = bincenters{2};

% Compute sensory noise std per trial for vision
sigmas_vis = sigma_fun(bincenters_vis,sigmazero_vis,w_vis);
if isscalar(sigmas_vis); sigmas_vis = sigmas_vis*ones(numel(bincenters_vis),1); end

% Compute sensory noise std per trial for vestibular
sigmas_vest = sigma_fun(bincenters_vest,sigmazero_vest,w_vest);
if isscalar(sigmas_vest); sigmas_vest = sigmas_vest*ones(numel(bincenters_vest),1); end

% Measurements
xrange_vis = zeros(1, grid_size, 1);
xrange_vest = zeros(1, 1, grid_size);
xrange_vis(1, :, 1) = alpha_rescaling_vis*linspace(max(min(bincenters_vis-MAXSD*sigmas_vis),-MAXRNG_XMEAS), min(max(bincenters_vis+MAXSD*sigmas_vis), MAXRNG_XMEAS), grid_size);
xrange_vest (1, 1, :) = alpha_rescaling_vest*linspace(max(min(bincenters_vest-MAXSD*sigmas_vest),-MAXRNG_XMEAS), min(max(bincenters_vest+MAXSD*sigmas_vest), MAXRNG_XMEAS), grid_size);
dx_vis = xrange_vis(1, 2, 1) - xrange_vis(1, 1, 1);
dx_vest = xrange_vest(1, 1, 2) - xrange_vest(1, 1, 1);

% Wrap large noisy measurement around circle?
wraparound_flag = MAXRNG_XMEAS >= 180 && ...
        ( min(bincenters_vis-MAXSD*sigmas_vis) <= -180 || max(bincenters_vis+MAXSD*sigmas_vis) >= 180 || ...
        min(bincenters_vest-MAXSD*sigmas_vest) <= -180 || max(bincenters_vest+MAXSD*sigmas_vest) >= 180);

srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';
ds = diff(srange(1:2));

if nargout > 1 % Save variables for debug or data generation
    extras.xrange_vis = xrange_vis;
    extras.xrange_vest = xrange_vest;
    extras.srange = srange;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute likelihood; variable range depends on type

likerange_vis = srange; 
likerange_vest = srange; 

sigmalikezero_vis = sigmazero_vis;
sigmalikezero_vest = sigmazero_vest;
wlike_vis = w_vis;
wlike_vest = w_vest;

% Compute sensory likelihood std for vision
sigmasprime_vis = sigma_fun(likerange_vis,sigmalikezero_vis,wlike_vis);

% Compute sensory likelihood std for vestibular
sigmasprime_vest = sigma_fun(likerange_vest,sigmalikezero_vest,wlike_vest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute likelihood for non-Gaussian likelihoods

if wraparound_flag
    like_vis = bsxfun_normpdf(xrange_vis,srange,sigmasprime_vis) + ...
        bsxfun_normpdf(xrange_vis,srange + 360,sigmasprime_vis) + ...
        bsxfun_normpdf(xrange_vis,srange - 360,sigmasprime_vis);
    like_vest = bsxfun_normpdf(xrange_vest,srange,sigmasprime_vest) + ...
        bsxfun_normpdf(xrange_vest,srange + 360,sigmasprime_vest) + ...
        bsxfun_normpdf(xrange_vest,srange - 360,sigmasprime_vest);
else
    like_vis = bsxfun_normpdf(xrange_vis,srange,sigmasprime_vis);
    like_vest = bsxfun_normpdf(xrange_vest,srange,sigmasprime_vest);
end

% Compute UNCORRELATED prior, p(s)
priorpdf1d = bsxfun_normpdf(srange,prior_mu,prior_sigma);
priorpdf1d = priorpdf1d/(qtrapz(priorpdf1d, 1)*ds); % Normalize prior

% Compute unnormalized posterior and rightward posterior (C = 2)
postpdf_c2 = bsxfun(@times, priorpdf1d, like_vest);

postright_c2 = [];
if do_estimation
    if priorc1 < 1
        postright_c2(1,:,:) = VestBMS_PostRight(postpdf_c2);
    else
        postright_c2 = 0;
    end
end

likec1 = [];

% Compute unnormalized posterior and rightward posterior (C = 1)
if priorc1 > 0
    if do_estimation
        [postright_c1(1,:,:),likec1(1,:,:)] = VestBMS_c1postandlikec1qtrapz(postpdf_c2, like_vis);
        likec1 = likec1*ds + realmin;
    else
        postright_c1 = [];
    end
else
    postright_c1 = 0;
end

if nargout > 1 
    extras.postright_c1 = postright_c1; 
    extras.postright_c2 = postright_c2;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute causal inference weights

% Compute marginal likelihood, p(x_vis, x_vest|C)

% CASE C=2, Independent likelihoods
likec2_vis = bsxfun(@times, priorpdf1d, like_vis);
likec2 = (bsxfun(@times, qtrapz(likec2_vis, 1)*ds, qtrapz(postpdf_c2, 1)*ds)) + realmin;

% postpdfc1 = bsxfun(@rdivide, postpdfc1, likec1);

% CASE C=1, Likelihoods are not independent
if isempty(likec1)
    likec1(1,:,:) = VestBMS_likec1qtrapz(postpdf_c2,like_vis)*ds + realmin;
end

% Compute weight for cue fusion
% Model weight is Bayesian posterior p(C=1|x_1, x_2)
w1 = likec1*priorc1./(likec1*priorc1 + likec2*(1-priorc1));

% NaNs can emerge as 0/0 - assume that the response becomes random
w1(isnan(w1)) = 0.5;

% Bisensory estimation
if do_estimation
    % Compute posterior probability of rightward motion    
    postright = bsxfun(@plus, bsxfun(@times, w1, postright_c1), bsxfun(@times, 1-w1, postright_c2));

    % Probability of rightward response            
    prright = 1./(1 + ((1-postright)./postright).^beta_softmax);
    prright(isnan(prright)) = 0.5;
end

% Bisensory unity judgement
if do_unity
    if beta_softmax == 1
        w1_unity = w1;
    elseif beta_softmax == Inf
        w1_unity = zeros(size(w1));
        w1_unity(w1 > 0.5) = 1;
        w1_unity(w1 == 0.5) = 0.5;
    else
        w1_unity = 1./(1 + ((1-w1)./w1).^beta_softmax);
    end
end

if nargout > 1 % Save variables for debug or data generation
    extras.w1 = w1;
    % extras.postpdfc1 = postpdfc1;
end

%----------------------------------------------------------------------
% Plot decision rule for unity judgments (function of x_vest, x_vis)
if DEBUG && do_unity
    clf;
    surf(xrange_vis(:), xrange_vest(:), squeeze(w1_unity),'LineStyle','none'); hold on;
    xlabel('$x_\mathrm{vis}$ (deg)','Interpreter','LaTeX','FontSize',14); 
    ylabel('$x_\mathrm{vest}$ (deg)','Interpreter','LaTeX','FontSize',14);
    axis square;
    xylims = [min(xrange_vest(1),xrange_vis(1)),max(xrange_vest(end),xrange_vis(end))];
    xlim(xylims);
    ylim(xylims);
    plot3(xylims,xylims,[200 200],'--k','LineWidth',1);
    view([0 90]);
    hold off;
    title(['$\sigma^0_\mathrm{vis} = ' num2str(sigmazero_vis,'%.2f') '$ deg, $w_\mathrm{vis} = ' num2str(abs(w_vis),'%.2f') '$' ...
        ', $\sigma^0_\mathrm{vest} = ' num2str(sigmazero_vest,'%.2f') '$ deg, $w_\mathrm{vest} = ' num2str(abs(w_vest),'%.2f') '$'], ...
        'FontSize',14,'Interpreter','LaTeX');
    pause;
end
%----------------------------------------------------------------------    

% Clean up memory
clear postright w1 postpdf_c2 postright_c1 postright_c2 ...
    likec1 likec2 likec2_vis postpdfc1 lratio;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Marginalize over noisy measurements

xpdf_vis = bsxfun_normpdf(xrange_vis, alpha_rescaling_vis*bincenters_vis,alpha_rescaling_vis*sigmas_vis);
if wraparound_flag
    xpdf_vis = xpdf_vis + ...
        bsxfun_normpdf(xrange_vis, alpha_rescaling_vis*bincenters_vis + 360,alpha_rescaling_vis*sigmas_vis) ...
        + bsxfun_normpdf(xrange_vis, alpha_rescaling_vis*bincenters_vis - 360,alpha_rescaling_vis*sigmas_vis);        
end
xpdf_vis = bsxfun(@rdivide, xpdf_vis, qtrapz(xpdf_vis, 2)); % Not multiplying by volume element

xpdf_vest = bsxfun_normpdf(xrange_vest, alpha_rescaling_vest*bincenters_vest,alpha_rescaling_vest*sigmas_vest);    
if wraparound_flag
    xpdf_vest = xpdf_vest + ...
        bsxfun_normpdf(xrange_vest, alpha_rescaling_vest*bincenters_vest + 360,alpha_rescaling_vest*sigmas_vest) ...
        + bsxfun_normpdf(xrange_vest, alpha_rescaling_vest*bincenters_vest - 360,alpha_rescaling_vest*sigmas_vest);        
end    
xpdf_vest = bsxfun(@rdivide, xpdf_vest, qtrapz(xpdf_vest, 3));  % Not multiplying by volume element

if do_estimation
    prmat = zeros(numel(bincenters_vest), 2);
    prmat(:,2) = VestBMS_finalqtrapz(xpdf_vis,xpdf_vest,prright);   % Not multiplying by volume element (xpdfs did not)
    prmat(:,1) = 1 - prmat(:,2);
else
    prmat = [];
end

if do_unity
    prmat_unity = zeros(numel(bincenters_vest), 2);
    prmat_unity(:,1) = VestBMS_finalqtrapz(xpdf_vis,xpdf_vest,w1_unity);    % Not multiplying by volume element (xpdfs did not)
    prmat_unity(:,2) = 1 - prmat_unity(:,1);
else
    prmat_unity = [];
end

% Fix probabilities
prmat = min(max(prmat,0),1);
prmat_unity = min(max(prmat_unity,0),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalize log likelihood

prmat = lambda/2 + (1-lambda)*prmat;
prmat_unity = lambda/2 + (1-lambda)*prmat_unity;

[ll,extras] = finalize(prmat,prmat_unity,counts,FIXEDLAPSEPDF,nargout > 1,sum_flag);
varargout{1} = ll;
if nargout > 1; varargout{2} = extras; end

end

%--------------------------------------------------------------------------
function [ll,extras] = finalize(prmat,prmat_unity,counts,epsilon,extras_flag,sum_flag)
%FINALIZE Finalize log likelihood

prmat = 0.5*epsilon + (1-epsilon)*prmat;
prmat_unity = 0.5*epsilon + (1-epsilon)*prmat_unity;

if extras_flag
    extras.responsepdf = prmat;
    extras.responsepdf_unity = prmat_unity;
else
    extras = [];
end

prmat = [prmat(:); prmat_unity(:)];
xx = [counts{1}(:); counts{2}(:); counts{3}(:)];

if sum_flag
    ll = sum(xx.*log(prmat));
else
    ll = loglikmat2vec(log(prmat),xx);
end

end
%--------------------------------------------------------------------------
function postright = VestBMS_PostRight(postpdf)
%VESTBMS_POSTRIGHT Posterior probability of perceiving rightward motion.

n = size(postpdf,1);
pmin = realmin*n/2;

idx0deg = (n+1)/2; % Index of 0 deg

% Use MEX files
% postleft = qtrapzc(postpdf,1,[1;idx0deg]) + pmin;
% posttemp = qtrapzc(postpdf,1,[idx0deg;size(postpdf,1)]) + pmin;
% postright(1,:,:) = posttemp./(posttemp + postleft);
    
    %t1 = postright;
    %clear postright;

postleft = qtrapz(postpdf(1:idx0deg,:,:),1) + pmin;
postright = qtrapz(postpdf(idx0deg:end,:,:),1) + pmin;
postright(1,:,:) = postright./(postright + postleft);
    
    %sum(abs(t1(:) - postright(:)))

end
%--------------------------------------------------------------------------
function [postright_c1,likec1] = VestBMS_c1postandlikec1qtrapz(postpdf_c2,like_vis)
%VESTBMS_C1POSTANDLIKEC1QTRAPZ Multiple computations for C=1 (uncorrelated)
%
% ================ INPUT VARIABLES ====================
% POSTPDF_C2: p(s) * p(x_vest|s). [S-by-1-by-K] (double)
% LIKE_VIS: p(x_vis|s). [S-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% POSTRIGHT_C1: p(right|x_vis,x_vest,C=1). [K-by-K] (double)
% LIKEC1: p(x_vis,x_vest|C=1). [K-by-K] (double)

postpdf_c1 = bsxfun(@times,postpdf_c2,like_vis);
postright_c1(:,:) = VestBMS_PostRight(postpdf_c1);
likec1(:,:) = qtrapz(postpdf_c1,1);
end
%--------------------------------------------------------------------------
function likec1 = VestBMS_likec1qtrapz(postpdf_c2,like_vis)
%VESTBMS_LIKEC1QTRAPZ Compute p(x_vis,x_vest|C=1) for uncorrelated prior
%
% ================ INPUT VARIABLES ====================
% POSTPDF_C2: p(s) * p(x_vest|s). [S-by-1-by-K] (double)
% LIKE_VIS: p(x_vis|s). [S-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% LIKEC1: p(x_vis,x_vest|C=1). [K-by-K] (double)

likec1(:,:) = qtrapz(bsxfun(@times,postpdf_c2,like_vis),1);
end
%--------------------------------------------------------------------------
function prmat = VestBMS_finalqtrapz(xpdf_vis,xpdf_vest,R)
%VESTBMS_FINALQTRAPZ Marginalize response probability over x_vis and x_vest
%
% ================ INPUT VARIABLES ====================
% XPDF_VIS: p(x_vis|s). [S-by-K] (double)
% XPDF_VEST: p(x_vest|s). [S-by-1-by-K] (double)
% R: p(response|x_vis,x_vest). [1-by-K-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% RPDF: p(r|s). [S-by-1] (double)

prmat = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_vis, xpdf_vest), R), 2), 3);
end
%--------------------------------------------------------------------------
function z = qtrapz(y,dim)
%QTRAPZ  Quick trapezoidal numerical integration.
%   Z = QTRAPZ(Y) computes an approximation of the integral of Y via
%   the trapezoidal method (with unit spacing).  To compute the integral
%   for spacing different from one, multiply Z by the spacing increment.
%
%   For vectors, QTRAPZ(Y) is the integral of Y. For matrices, QTRAPZ(Y)
%   is a row vector with the integral over each column. For N-D
%   arrays, QTRAPZ(Y) works across the first non-singleton dimension.
%
%   Z = QTRAPZ(Y,DIM) integrates across dimension DIM of Y. The length of X 
%   must be the same as size(Y,DIM).
%
%   QTRAPZ is up to 3-4 times faster than TRAPZ for large arrays.
%
%   See also TRAPZ.

% Luigi Acerbi <luigi.acerbi@nyu.edu>
% Version 1.0. Release date: Jul/20/2015.

% By default integrate along the first non-singleton dimension
if nargin < 2; dim = find(size(y)~=1,1); end    

% Behaves as sum on empty array
if isempty(y); z = sum(y,dim); return; end

% Compute dimensions of input matrix    
if isvector(y); n = 1; else n = ndims(y); end

switch n
    case {1,2}      % 1-D or 2-D array
        switch dim
            case 1
                z = sum(y,1) - 0.5*(y(1,:) + y(end,:));
            case 2
                z = sum(y,2) - 0.5*(y(:,1) + y(:,end));
            otherwise
                error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
        end

    case 3      % 3-D array
        switch dim
            case 1
                z = sum(y,1) - 0.5*(y(1,:,:) + y(end,:,:));
            case 2
                z = sum(y,2) - 0.5*(y(:,1,:) + y(:,end,:));
            case 3
                z = sum(y,3) - 0.5*(y(:,:,1) + y(:,:,end));
            otherwise
                error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
        end                

    case 4      % 4-D array
        switch dim
            case 1
                z = sum(y,1) - 0.5*(y(1,:,:,:) + y(end,:,:,:));
            case 2
                z = sum(y,2) - 0.5*(y(:,1,:,:) + y(:,end,:,:));
            case 3
                z = sum(y,3) - 0.5*(y(:,:,1,:) + y(:,:,end,:));
            case 4
                z = sum(y,4) - 0.5*(y(:,:,:,1) + y(:,:,:,end));
            otherwise
                error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
        end                

    otherwise   % 5-D array or more
        for iDim = 1:n; index{iDim} = 1:size(y,iDim); end
        index1 = index;     index1{dim} = 1;
        indexend = index;   indexend{dim} = size(y,dim);
        try
            z = sum(y,dim) - 0.5*(y(index1{:}) + y(indexend{:}));
        catch
            error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');            
        end
end

end
%--------------------------------------------------------------------------
function y = bsxfun_normpdf(x,mu,sigma)
%BSXFUN_NORMPDF Vectorized normal probability density function (pdf).
%   Y = BSXFUN_NORMPDF(X,MU,SIGMA) returns the pdf of the normal 
%   distribution with mean MU and standard deviation SIGMA, evaluated at 
%   the values in X. Dimensions of X, MU, and SIGMA must either match, or 
%   be equal to one. Computation of the pdf is performed with singleton
%   expansion enabled via BSXFUN. The size of Y is the size of the input 
%   arguments (expanded to non-singleton dimensions).
%
%   All elements of SIGMA are assumed to be non-negative (no checks).
%
%   See also BSXFUN, BSXFUN_NORMCDF, NORMPDF.

%   Author: Luigi Acerbi
%   Release date: 15/07/2015

if nargin<3
    error('bmp:bsxfun_normpdf:TooFewInputs','Input argument X, MU or SIGMA are undefined.');
end

try
    if isscalar(mu)
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, x - mu, sigma).^2), sigma)/sqrt(2*pi);
    elseif isscalar(sigma)
        y = exp(-0.5*(bsxfun(@minus, x, mu)/sigma).^2)/(sigma*sqrt(2*pi));
    else
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma).^2), sigma)/sqrt(2*pi);
    end
catch
    error('bmp:bsxfun_normpdf:InputSizeMismatch',...
          'Non-singleton dimensions must match in size.');
end

end
%--------------------------------------------------------------------------
function [bimodal_counts,bincenters_bim] = data2counts(data,bincenters)
%DATA2COUNTS Bin bimodal data
% Prepare stimulus/response bin counts

num_noise = 3; % Number of noise conditions

bincenters_vis = repmat(bincenters, [length(bincenters), 1]);
bincenters_vest = repmat(bincenters', [length(bincenters), 1]);
bincenters_vest = bincenters_vest(:);
binmeshLength = size(bincenters_vis,1);

xmat_tot = zeros(binmeshLength,1);

for iNoise = 1:num_noise
    bimodal_counts{iNoise}{1} = [];
    for iType = 2:3
        if iType == 3
            responses = [1,2]; % Unity
        else
            responses = [1,-1]; % Estimation
        end
        
        data_now = data(data(:,3) == iType & data(:,7) == iNoise, :);
        
        xmat = zeros(binmeshLength,numel(responses));
        for iBin = 1:binmeshLength
            for jBin = 1:numel(responses)
                xmat(iBin,jBin) = sum(data_now(:, 5) == bincenters_vis(iBin) & ...
                data_now(:, 4) == bincenters_vest(iBin) & ...
                    data_now(:, 6) == responses(jBin));
            end
        end
        xmat_tot = xmat_tot + sum(xmat,2);
        bimodal_counts{iNoise}{iType} = xmat;
    end
end

idx_bin = xmat_tot > 0;

bincenters_bim{1} = bincenters_vis(idx_bin);
bincenters_bim{2} = bincenters_vest(idx_bin);

for iNoise = 1:num_noise
    for iType = 2:3
        bimodal_counts{iNoise}{iType} = bimodal_counts{iNoise}{iType}(idx_bin,:);
    end
end


end

