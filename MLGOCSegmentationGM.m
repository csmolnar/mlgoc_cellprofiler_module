function finalPhi = MLGOCSegmentationGM(PriorPhasefieldParameters, ExtendedImage, DataParameters, Kappa, ExtendedInitialPhi, OptimizationParameters)
%MLGOCSEGMENTATION Gradient descent optimization of multi layer 'gas of
%   circles model and additive data model.
%   MLGOCSegmentationGM(PriorPhasefieldParameters, ExtendedImage, DataParameters, Kappa, ExtendedInitialPhi, OptimizationParameters)
%   returns a multi layered phase field as a local minimum of the energy
%   the multi layer 'gas of circles phase field model combined with
%   additive data image model after gradient descent optimization. The
%   optimizer start the phase field at ExtendedInitialPhi.

[hExtended, wExtended] = size(ExtendedImage);

LayerNumber = length(PriorPhasefieldParameters);

maxd = int32( max( [PriorPhasefieldParameters.d] ) );

hImage = hExtended - 2*maxd;
wImage = wExtended - 2*maxd;

threshold = zeros(LayerNumber,1);
for ll=1:LayerNumber
    threshold(ll) = PriorPhasefieldParameters(ll).alpha/PriorPhasefieldParameters(ll).lambda;
end

fPhi = mlEvolution(ExtendedInitialPhi, Kappa, OptimizationParameters.tolerance, OptimizationParameters.maxIts, PriorPhasefieldParameters, DataParameters, OptimizationParameters.saveFreq, ExtendedImage);

finalPhi = fPhi(maxd+1:maxd+hImage,maxd+1:maxd+wImage,:);

end

function newPhi = mlEvolution(initPhi, kappa, tolerance, maxIts, parameters, dataParameters, showFreq, extImage)
%MLEVOLUTION Gradient descend method implementation

% useTanh=1 makes the computations more stable
useTanh = 1;

% Fourier transform of the kernel
linearOp = computeLinearPart(initPhi, parameters);

% imageLinearOp = computeImageLinearPart(image, parameters);

%% Initialize loop
% layerNum = size(initPhi, 3);

oldPhi = initPhi;

converged = 0;

numIts = 0;

% deltaPhi = 0;

% Ep = zeros(maxIts,layerNum);
% Ed = zeros(maxIts,1);
% 
% EdDenomSum = zeros(maxIts,1);


newPhi = initPhi;

% meanFuncDer = 0;
% meanDeltaPhi = 0;
% maxFuncDer = 0;
% deltaT = 0;

% threshold = zeros(layerNum, 1);

% for i=1:layerNum
%     threshold(i) = parameters(i).alpha/parameters(i).lambda;
% end

%totenergy = zeros(1, 10000);
% m = 0;

means = zeros(maxIts,1);
maxs = zeros(maxIts,1);

%% Optimalization
while (~converged)
    
%     if showFreq>0
%         if (mod(numIts, showFreq) == 0)
% %             toc
%             
%             name = ['out' num2str(numIts, '%07d')];
%             
% %             parfor l=1:layerNum
%             for l=1:layerNum
%                 actPhi(:,:,l) = oldPhi(parameters(1).d+1:parameters(1).d+Nx,parameters(1).d+1:parameters(1).d+Ny,l);
%             end
%             
%             coi = createContourImage(['input\' parameters(1).imageName], actPhi, threshold);
%             coiName = ['iteration_' fileName secondPart name];
%             coiName(coiName=='.') = '_';
%             coiName = [coiName '_coi.bmp']
%             imwrite(coi, coiName);
%             
%         end
%         
%     end
    
    [funcDer, overlapDer] = mlEvolveStep(oldPhi, linearOp, parameters, dataParameters, extImage, kappa, useTanh);
    
    funcDer = funcDer + overlapDer;
    
    maxFuncDer = max(abs(funcDer(:)));
    
    maxs(numIts+1) = maxFuncDer;
    
%     deltaT = 0.025;
    
    deltaT = 1/(10*maxFuncDer);
    
    deltaPhi = -deltaT*funcDer;
    
    newPhi = oldPhi + deltaPhi;
    
    meanFuncDer = mean( abs( funcDer(:) ) );
    means(numIts+1) = meanFuncDer;
   
%     if mod(numIts, 100) == 0
%        fprintf('Iteration: %d\r', numIts);
%         numIts
%         meanFuncDer
%         tolerance
%         if layerNum==1 && 0
%             edgeImage = edge(newPhi>(parameters(1).alpha/parameters(1).lambda));
%             imageEdge = extImage;
%             imageEdge(edgeImage) = 1;
%             figure(1); imagesc(imageEdge); colorbar; colormap(gray); title([num2str(numIts) '. iteration']);
%         end
%     end
    
    
%     [E, ~, ~, ~] = mlEnergy(oldPhi, parameters);
%     Ep(numIts) = sum(E);

%% calculate energies 
%     for ll=1:layerNum
%         Ep(numIts, ll) = energy(oldPhi(:,:,ll), parameters(ll));
%     end
%     
%     if useTanh ~= 1
%  
%         Ed(numIts) = dataEnergyProstate(oldPhi, image, parameters);
% 
%         EdDenomSum(numIts) = computeDenomDataSum(oldPhi, parameters);
% 
%     else
%     
%         Ed(numIts) = dataEnergyProstateTanh(oldPhi, image, parameters);
%         
%         EdDenomSum(numIts) = computeDenomDataSumTanh(oldPhi, parameters);
%     
%     end
    
    if mod(numIts, int16(maxIts/10))==0
        fprintf('Iteration %6d (%3d%%)\n', numIts, int16(100.0*numIts/maxIts));
    end
    
%% change phase field

    oldPhi = newPhi;
    
    numIts = numIts + 1;
    
    if ((meanFuncDer < tolerance) || (numIts >= maxIts))
        converged = 1;
        fprintf('Iteration %6d (%3d%%)\n', maxIts, 100);
    end
    
end

end

function linearOp = computeLinearPart(initPhi, parameters)
%COMPUTELINEARPART 

[h, w, layerNum] = size(initPhi);

linearOp = cell(layerNum,1);

for i = 1:layerNum
    
    D = parameters(i).D;
    
    lambda = parameters(i).lambda;
    
    beta = parameters(i).beta;
    
    discrete = parameters(i).discrete;
    
    K2 = computeNegLaplacian([h w], discrete);
    
    interactionOperator = computeInteractionOperator(K2, parameters(i));
    
    linearOp{i} = K2 .* (D - beta*interactionOperator) - lambda;
    
    %     Code added to facilitate filtering in real space not Fourier space.
    
%     if (parameters(i).boundarycorr == 1)
%         
%         linearOp{i} = real(fftshift(ifft2(ifftshift(linearOp{i}))));
%         
%         % Crop linearOp to support.
%         linearOp{i} = linearOp{i}([Nx/2 + 1 - 2*d:Nx/2 + 1 + 2*d], [Ny/2 + 1 - 2*d:Ny/2 + 1 + 2*d]);
%         
%     end
    
end

end

function K2 = computeNegLaplacian(N, discrete)

Ny = N(1);

Nx = N(2);


kx = -pi:(2*pi/Nx):pi*(Nx - 1)/Nx;

ky = -pi:(2*pi/Ny):pi*(Ny - 1)/Ny;


[Kx, Ky] = meshgrid(kx, ky);

K2d = 2*(1 - cos(Kx) + 1 - cos(Ky));

K2c = Kx .* Kx + Ky .* Ky;

if (discrete)
    
    K2 = K2d;
    
else
    
    K2 = K2c;
    
end

end

function interactionOperator = computeInteractionOperator(K2, parameters)

if parameters.marie~=1
% if ~strcmp(parameters.funName, 'Marie')
    
    interactionOperator = 1 ./(1 + K2);
    
else
    
    N = size(K2);
    
    Nx = N(2);
    
    Ny = N(1);
    
    d = parameters.d;
    
    epsilon = parameters.epsilon;
    
    
    x = -Nx/2:1:Nx/2 - 1;
    
    y = -Ny/2:1:Ny/2 - 1;
    
    
    [X, Y] = meshgrid(x, y);
    
    R = sqrt(X .* X + Y .* Y);
    
    innerSpaceIntInterOp = (R <= d - epsilon);
    
    centredR = (R - d) / epsilon;
    
    centreSpaceIntInterOp = (1/2) * (1 - centredR - (1/pi) * sin(pi * centredR)) .* (R > d - epsilon) .* (R < d + epsilon);
    
    spaceIntInterOp = innerSpaceIntInterOp + centreSpaceIntInterOp;
    
    interactionOperator = real(fftshift(fft2(ifftshift(spaceIntInterOp))));
    
end

end

function nonlinearPart = computeNonlinearPart(phi, lambda, alpha)

phi2 = phi.*phi;

phi3 = phi.*phi2;

nonlinearPart = lambda*phi3 - alpha*phi2;

end

function [funcDer, overlapDer] = mlEvolveStep(oldPhi, linearOp, parameters, dataParameters, image, kappa, useTanh)

% Note that the Phi variables are in Fourier space.
% linearPart is (a multiple of) the Fourier transform of Laplacian, passed in to make things quicker.
% The zero of Fourier space is at (need to check how FFT works in scilab).

[~,~,layerNum] = size(oldPhi);

kappas = zeros(layerNum,1) + kappa;

funcDer = zeros(size(oldPhi));
overlapDer = zeros(size(oldPhi));

sumPhi = sum(oldPhi,3);

% compute the product of the layers and then the use it with data
% imagePart = computeImagePart(oldPhi, image, parameters);

if useTanh ~= 1
    
    % simple sum
    imagePart = computeImagePartAdditive(oldPhi, image, dataParameters);
    
else
    
    % sum_i {(tanh(phi^i)+1)/2}
    imagePart = computeImagePartAdditiveTanh(oldPhi, image, dataParameters);
    
end

for i=1:layerNum
    
    % Linear factor
    
    % If boundary correction is needed because of the Fourier transform
    
    %     if (parameters(i).boundarycorr == 0) % works in Fourier space
    
    hatOldPhi = fftshift(fft2(oldPhi(:,:,i)));
    if iscell(linearOp)
        opOldPhi = linearOp{i}.*hatOldPhi;
    else
        opOldPhi = linearOp.*hatOldPhi;
    end
    
    linearPart = real(ifft2(ifftshift(opOldPhi)));
    %     else
    %         %         works in image space
    %         %         linearPart = imfilter(oldPhi, linearOp, 'replicate');
    %     end
    
    % Compute the nonlinear contribution, which can be zero.
    
    nonlinearPart = computeNonlinearPart(oldPhi(:,:,i), parameters(i).lambda, parameters(i).alpha);
    
    funcDer(:,:,i) = linearPart + parameters(i).alpha + nonlinearPart + imagePart(:,:,i);
    
    %     Add functional derivative of overlapping energy
    if layerNum>1
        overlapDer(:,:,i) = kappas(i)/2 * (sumPhi - oldPhi(:,:,i) + (layerNum-1) );
    end
    
end

end

function imageLinearPart = computeImagePartAdditive(oldPhi, image, dataParameters)
% computeImagePartProstate computes the functional derivative of data term of
% energy

[~, ~, layerNum] = size(oldPhi);

phiPlus = sum(oldPhi, 3)/2 + layerNum/2;
% phiPlus = sum((oldPhi+1)/2, 3);

phiPlus2 = phiPlus.^2;

[Nx, Ny] = size(image);

imageLinearPart = zeros(Nx, Ny, layerNum);

gamma2 = dataParameters.gamma2;
    
%     gamma1 = parameters(i).gamma1;
    
    muin = dataParameters.muin;
    
    muout = dataParameters.muout;
    
%     weightsin = parameters(i).win;
    
    sigmain = dataParameters.sigmain;
    
    sigmaout = dataParameters.sigmaout;
    
%     weightout = parameters(i).wout;

for i=1:1
    
%     intensityPart = zeros(Nx, Ny);
    
    deltaMu = muin-muout;
    deltaMu2 = deltaMu^2;
    deltaSigma2 = sigmain^2-sigmaout^2;
    sigmaout2 = sigmaout^2;
    
    intensityPart = (deltaMu2*deltaSigma2) * phiPlus2 ...
                     + (2*sigmaout2*deltaMu2) * phiPlus ...
                     - (2*sigmaout2*deltaMu) * (image-muout) ...
                     - deltaSigma2 * (image-muout).^2;
    
    intensityPart = intensityPart ./ ( sigmaout2 + (deltaSigma2) * phiPlus ).^2;
    
    imageLinearPart(:,:,i) = (gamma2 /4 ) * intensityPart;
    
end

% The functional derivative is the same for all layers:

if layerNum>1
    for i=2:layerNum
        imageLinearPart(:,:,i) = imageLinearPart(:,:,1);
    end
end

end

function imageLinearPart = computeImagePartAdditiveTanh(oldPhi, image, dataParameters)
% computeImagePartProstate computes the functional derivative of data term of
% energy

[~, ~, layerNum] = size(oldPhi);

tanhPhi = tanh(oldPhi);
tildePhiPlus = sum(tanhPhi, 3) + layerNum/2;
tildePhiPlus2 = tildePhiPlus.^2;

sech2Phi = sech(oldPhi).^2;

[Nx, Ny] = size(image);

imageLinearPart = zeros(Nx, Ny, layerNum);

gamma2 = dataParameters.gamma2;

%     gamma1 = parameters(i).gamma1;

muin = dataParameters.muin;

muout = dataParameters.muout;

%     weightsin = parameters(i).win;

sigmain = dataParameters.sigmain;

sigmaout = dataParameters.sigmaout;

%     weightout = parameters(i).wout;

deltaMu = muin-muout;
deltaMu2 = deltaMu^2;
deltaSigma2 = sigmain^2-sigmaout^2;
sigmaout2 = sigmaout^2;

intensityPartNominator = (deltaMu2*deltaSigma2) * tildePhiPlus2 ...
    + (2*sigmaout2*deltaMu2) * tildePhiPlus ...
    - (2*sigmaout2*deltaMu) * (image-muout) ...
    - deltaSigma2 * (image-muout).^2;
intensityPartDenominator = ( sigmaout2 + (deltaSigma2) * tildePhiPlus2 ).^2;
intensityPart = ( gamma2/4 ) * intensityPartNominator ./ intensityPartDenominator;

for i=1:layerNum
    imageLinearPart(:,:,i) =  intensityPart .* sech2Phi(:,:,i); 
end

end