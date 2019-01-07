function [psf, OTF, csfLambda] = calcPsfMicroscope(psfSz, voxSz, lambdas, varargin)
% 2017.04.27_Xiang & Kuan modified
% psf: conventional psf
% restrict lambdas to odd numbers, such as lambdas =
% linspace(515*(10^(-3)),525*(10^(-3)),11) for multiple wavelength; or
% lambdas = 520*(10^(-3)) for single wavelength
% wLambdas =  the weight for each lambda
% different lambdas may have different weight
% total number of psf slice is: 
% % psf outside microscope (2*N1+1), verify psf at the border is very weak, otherwise increase N1
% voxSz:    the psf voxSz in object space, nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% after this run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MicroscopeMFM3DPSF.m on the
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% left by changing voxSz(2) = 0.05;
%% MFM DOE load
% a1 + a2 = circlular support
%% Setup Parameters
%% Parse positional Parameters
if ~exist('psfSz', 'var') || isempty(psfSz)    
    psfSz = [101, 101, 121];    
end
if ~exist('lambdas', 'var') || isempty(lambdas)    
    lambdas = 520;    
    %lambdas = linspace(515*(10^(-3)),525*(10^(-3)),11);
end
NLambda = length(lambdas);
%% Parse name-value parameters
p = inputParser;  % CaseSensitive=0, PartialMatching=1
p.addParameter('verbose',   true); % true/false, whether show result
p.addParameter('wLambdas',  ones(1,NLambda),    @(x) isempty(x) || (isnumeric(x) && isvector(x))); % 1D numeric vector
% Microcope Parameters
p.addParameter('M',         60); % maginification size; 20x/0.5 or 40x/0.95 
p.addParameter('NA',        1.27); % 
p.addParameter('f_tl',      200*1e6); % um, tube lens 200mm 
p.parse(varargin{:});
cellfun(@(x) evalin('caller', sprintf('%s=p.Results.%s;', x, x)), fields(p.Results)); % e.g. eval('fSz = p.Results.fSz'); Or use in.fSz after set in =p.Results;
% compute dependent parameters
if  ~exist('voxSz', 'var') || isempty(voxSz)
    voxSz(1:2) = 13*1e3/(400/200)/M; %camera_pixel_size/(f2/f1*M), camera_pixel_size=13um, f2=400mm, f1=200mm; Pixel size after Microscope ,why much less than  diff_limit_x=central_lambda/2/NA= 0.205um    
    voxSz(3) = 50; % 50nm, much less thant diff_limit_z = 2*central_lambda/(NA)^2 = 0.645um;
end
N1 = floor(psfSz(1)/2);
f_obj = f_tl/M; % objective lens
n = NA/0.95; % must be > NA, water refractive index ??? questionable? microscope wa


%% Voxel size in the object space, also compare with diffraction limits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Voxel size after microscope: code from dorgon nature LFM code
x1space = voxSz(1)*M*(-N1:1:N1); 
x2space = voxSz(2)*M*(-N1:1:N1);
psf = zeros(psfSz);
p1 = 0;  % y of psf(u, v; x, y, z);
p2 = 0;  % x of psf(u, v; x, y, z);
%pad_num = (size(a2, 1) - (2*N1+1))/2; % padd to same size as DOE

psfLambdaZ = zeros([psfSz(1:2), NLambda]);
csfLambda = zeros([psfSz, NLambda]);
%% 3D PSF generation
for nZ = 1:psfSz(3)
    p3 = (nZ-ceil((psfSz(3)+1)/2)) * voxSz(3);  % z of psf(u, v; x, y, z); from -50:50 etc; %zCenter = ceil((psfSz(3)+1)/2)) = 51
    for nLambda = 1:NLambda
        lambda = lambdas(nLambda); % um
       %% generate the field in native image plane u_i
        tic;
        csfLambda(:, :, nZ, nLambda) = calcPSF(p1, p2, p3, f_obj, NA, x1space, x2space, lambda,  M, n);
        psfLambdaZ(:, :, nLambda) = csfLambda(:,:,nZ,nLambda).*conj(csfLambda(:,:,nZ,nLambda));
        toc;
    end
    psf(:,:,nZ) = sum(psfLambdaZ,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psf = psf./(max(psf(:)));
if verbose
    figure;show3dPointsWithProjections(psf, voxSz, 0.001);title('3D microscope PSF')
end
OTF = (fftn(psf));
OTF = fftshift(abs(OTF));
OTF = OTF./max(OTF(:));
if verbose
    figure;show3dPointsWithProjections(OTF, 1./(voxSz.*psfSz), 0.01); title('3D microscope OTF')
end
end

