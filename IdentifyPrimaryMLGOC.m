function handles = IdentifyPrimaryMLGOC(handles)

% Help for the MLGOC Segmentation module:
% Category: Object Processing
%
% SHORT DESCRIPTION:
% Implementation of the multi-layered phase field 'Gas of Circles'
% variational model combined with an additive data model for e.g.
% fluorescent microscopy.
% *************************************************************************
%
% This module does the segmentation of touching or overlapping
% near-circular objects,  e.g. nuclei. The module can use the prelimitary
% results of IdentifyPrimAutomatic or DetectSpot module as initialization,
% but also can use initialization that does not require any initial object.
% The multi-layered 'gas of circles' shape model prefers near-circular
% objects with approximately same size. The corresponding image data model
% is an additive model assuming intensity of overlapping objects is
% proportional to the number of objects.
% 
% Molnar, C., Kato, Z., Jermyn, I. (2015). A New Model for the Segmentation
% of Multiple, Overlapping, Near-Circular Objects. DICTA 2015.
%
% Settings:
%
% Preferred radius of objects (in pixels):
% This parameters controls the interaction distance of 'gas of circles'
% model. The model prefers near-circular objects with the given radius.
%
% Number of layers: 
% Sets the number of layers of the multi-layer 'gas of circles' model.
% Automatic will set 4 layers in all cases of initialization. In case of
% Neutral and Squares initialization we can give the number of layers as an
% integer number. Running time is proportional to the number of layers.
%
% Mean intensity of a single object (muin - positive real number):
% Parameter of image model - measured intensity of a single object.
% 
% Variance of intensities over single objects (sigmain - positive real
% number):
% 
% Mean intensity of background (muout - positive real number):
% Parameter of image model - measured intensity of the background.
% 
% Variance of intensities of background (sigmaout - positive real
% number):
%
% Weight of data term (positive real number):
% Relative weight of image model compared to the geometric term. Lower
% values force the model more circular objects, higher values effect the
% model fits more to image data. NOTE: values close to zero result empty
% segmentation because of the 'inflection point' parameter setting of the '
% circles model.
%
% Type of initialization:
% * Seeds (manual) - This initialization is based on the previous module's
% identification of primary objects.
% * Circular Seeds (manual) - This initialization is based on the previous 
% module's identification of primary objects by placing circles with radius
% half of the preferred radius is set in the options.
% * Neutral - This option does not require any initial objects, but can
% result objects splitted into parts.
% * Squares - This option does not require any initial objects, but can
% result objects splitted into parts.
%
% Maximum number of iterations:
% Sets the maximum number of iterations of the gradient descent
% optimization method. Running time is proportional to the iteration
% number.
% Discard objects touching borders:
% 
% Discard objects touching the border of the image:
% You can choose to discard objects that touch the border of the image.
% This is useful in cases when you do not want to make measurements of
% objects that are not fully within the field of view (because, for
% example, the area or the shape would not be accurate).
% 
% Discard overlapping objects:
% You can choose to discard objects that have overlapping with another
% object. You can leave out cells if you don't want to make possible
% inaccurate measurement over areas belonging multiple objects. 
% 
% See also Identify primary Identify Secondary modules.
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by Csaba Molnar 2015.
%

% $Revision: 1017 $

%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
%drawnow

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What did you call the image you want to process?
%infotypeVAR01 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,1});
%inputtypeVAR01 = popupmenu

%textVAR02 = What did you call the initial seeds for correction?
%infotypeVAR02 = objectgroup
SeedName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = What do you want to call the objects identified by this module?
%defaultVAR03 = CorrectedNuclei
%infotypeVAR03 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,3});

%textVAR04 = Preferred radius of objects (in pixels)?
%defaultVAR04 = 20
Radius = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,4}));

%textVAR05 = Number of layers?
%defaultVAR05 = Automatic
LayerNumString = char(handles.Settings.VariableValues{CurrentModuleNum,5});

%textVAR06 = Mean intensity of a single object (muin - positive real number)?
%defaultVAR06 = 0.75
muin = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,6}));

%textVAR07 = Variance of intensities over single objects (sigmain - positive real
% number)?
%defaultVAR07 = 0.01
sigmain = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,7}));

%textVAR08 = Mean intensity of the background (muout - positive real number)?
%defaultVAR08 = 0.25
muout = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,8}));

%textVAR09 = Variance of intensities over background (sigmain - positive real
% number)?
%defaultVAR09 = 0.01
sigmaout = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,9}));

%textVAR10 = Weight of data term (positive real number)?
%defaultVAR10 = 1.0
DataWeight = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,10}));

%textVAR11 = Type of initialization?
%choiceVAR11 = Seeds (manual)
%choiceVAR11 = Circular Seeds (manual)
%choiceVAR11 = Neutral
%choiceVAR11 = Squares
Initialization = char(handles.Settings.VariableValues{CurrentModuleNum,11});
%inputtypeVAR11 = popupmenu

%textVAR12 = Maximum number of iterations?
%defaultVAR12 = 200
MaxIterations = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,12}));

%textVAR13 = Discard objects touching borders?
%choiceVAR13 = Yes
%choiceVAR13 = No
DiscardBorder = char(handles.Settings.VariableValues{CurrentModuleNum,13});
%inputtypeVAR13 = popupmenu

%textVAR14 = Discard overlapping objects if the Jaccard index of the overlapping area and the smaller object is over
%defaultVAR14 = 0.0
OverlapJIThreshold = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,14}));

%textVAR15 = What do you want to call the outlines of the corrected objects (optional)?
%defaultVAR15 = Do not save
%infotypeVAR15 = outlinegroup indep
SaveOutlines = char(handles.Settings.VariableValues{CurrentModuleNum,15});

%%%VariableRevisionNumber = 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRELIMINARY CALCULATIONS & FILE HANDLING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow

%%% Reads (opens) the image you want to analyze and assigns it to a variable,
%%% "OrigImage".
OrigImage = CPretrieveimage(handles,ImageName,ModuleName,'MustBeGray','CheckScale');
[ImageHeight, ImageWidth] = size(OrigImage);

%%% AlphaTilde a parameter for active contour model
AlphaTilde = 1;

%%% LamdbaTilde a scaling parameter for active contour model
LambdaTilde = 1;

%%% Rhatstar energy normalization parameter for inflection point model
%%% its values have to be between 0.69 and 0.78
%%% for positive alpha and beta values
% Rhatstar = 0.75;

%%% Rhatstar energy normalization parameter for 'gas of circle' model
Rhatstar = 1;

%%% ContourParameters contains the active contour, phase field and MRF
%%% parameters of the GOC inflection point model
% MRGOCIPMParameters = computeMRGOCIPMparameters(LambdaTilde, Radius, Rhatstar);

%%% ContourParameters contains the active contour, phase field and MRF
%%% parameters of the GOC model
MRGOCParameters = computeMRGOCparameters(AlphaTilde, LambdaTilde, Radius, Rhatstar);

%%% PriorPhasefieldParameters contains the phase field GOC parameters
PriorPhasefieldParameters = MRGOCParameters(2);
PriorPhasefieldParameters.discrete = 1;
PriorPhasefieldParameters.marie = 1;

%%% Get saved segmented objects from the Pipeline.
SegmentedGrayScaleObjects = CPretrieveimage(handles,['Segmented', SeedName],ModuleName,'DontCheckColor','DontCheckScale',size(OrigImage));

%%% We set a minimal intensity variance to provide numerical stability of
%%% the model
% SigmaThreshold = 0.00001;
% 
% muin = mean2(OrigImage(SegmentedGrayScaleObjects>0));
% sigmain = std2(OrigImage(SegmentedGrayScaleObjects>0));
% 
% if sigmain<SigmaThreshold
%     sigmain = SigmaThreshold;
% end
% 
% muout = mean2(OrigImage(SegmentedGrayScaleObjects==0));
% sigmaout = std2(OrigImage(SegmentedGrayScaleObjects==0));
% if sigmaout<SigmaThreshold
%     sigmaout = SigmaThreshold;
% end

MaskSegmentedObject = SegmentedGrayScaleObjects>0;
% SE = strel('disk', round(Radius*0.2), 0);
% MaskSegmentedObject = imerode(MaskSegmentedObject, SE);

%%% GradientWeight is the positive weight of image gradient term. It is not
%%% used in adaptive model (Molnar et al., DICTA 2015)
GradientWeight = 0.0;
DataParameters = struct('muin', muin, 'sigmain', sigmain, 'muout', muout, 'sigmaout',...
    sigmaout, 'win', 1, 'wout', 1, 'gamma1', GradientWeight, 'gamma2', DataWeight);

Maxd = int32(max([PriorPhasefieldParameters.d]));

%%% ExtendedImage is prepared to handle the side effects of built in
%%% Fourier transform.
ExtendedImage = extendImage(OrigImage, Maxd, DataParameters(1).muout);
[ExtendedImageHeight, ExtendedImageWidth] = size(ExtendedImage);

if ~isnumeric(LayerNumString)
    % sort objects to layers
    
    % InitialPhi = sortGrayscaleObjects2Layers(SegmentedGrayScaleObjects.*MaskSegmentedObject, Radius);
    [SegmentedGrayScaleObjects, n] = splitObjects(SegmentedGrayScaleObjects, max(SegmentedGrayScaleObjects(:)), Radius);
    InitialPhi = sortGrayscaleObjects2LayersColoring(SegmentedGrayScaleObjects.*MaskSegmentedObject, 4);
    
end

%%% Tolerance of gradient descent method
Tolerance  = 0.00001;

%%% Setting parameters of gradient descent optimization.
%%% maxIts field - maximum number of iterations
%%% saveFreq field - frequency of saving the actual phase field as a segmentation
%%% tolerance field - stop criteria for gradient descent method
OptimizationParameters = struct('maxIts', MaxIterations, 'saveFreq', -1, 'tolerance', Tolerance);

%%% Kappa is the weight of overlap penalty - not used with the adaptive
%%% model (Molnar et al., DICTA 2015)
Kappa = 0.0;

%%% Initialization can be Seeds (manual), Circular Seeds(Manual), Neutral or Squared
if strcmp(Initialization, 'Seeds (manual)')
%%% Seeds (manual) initialization is based on segmentation of previous
%%% module.

    LayerNumber = size(InitialPhi, 3);
    ExtendedInitialPhi = zeros(ExtendedImageHeight, ExtendedImageWidth, LayerNumber) - 1;
    ExtendedInitialPhi(Maxd+1:Maxd+ImageHeight, Maxd+1:Maxd+ImageWidth, :) = (InitialPhi*2-1)*0.1;
    
elseif strcmp(Initialization, 'Circular Seeds (manual)')
%%% Circular Seeds (manual) initialization is based on segmentation of previous
%%% module. For each segmented object we put a small circle of half radius
%%% of the preferred objects as initial object seeds.

    LayerNumber = size(InitialPhi, 3);
    InitialCircularPhi = createInitMLPhiCircles(InitialPhi, round(Radius/2));
    ExtendedInitialPhi = zeros(ExtendedImageHeight, ExtendedImageWidth, LayerNumber) - 1;
    ExtendedInitialPhi(Maxd+1:Maxd+ImageHeight, Maxd+1:Maxd+ImageWidth, :) = (InitialCircularPhi*2-1)*0.1;
    
elseif strcmp(Initialization, 'Neutral')
%%% Neutral initialization sets the initial phase field to a Gaussian white
%%% noise with mean of alpha_f/lambda_f and a tiny variance (no manual
%%% seeds are needed).

    if ~isempty(str2num(LayerNumString))
        if isnan(str2num(LayerNumString)) || (str2num(LayerNumString)<1)
            LayerNumber = 4;
        else
            LayerNumber = str2num(LayerNumString);
        end
    else
        LayerNumber = 4;
    end
    ExtendedInitialPhi = randn(ExtendedImageHeight,ExtendedImageWidth,LayerNumber) + ones(ExtendedImageHeight,ExtendedImageWidth,LayerNumber)*PriorPhasefieldParameters(1).alpha/PriorPhasefieldParameters(1).lambda;
    
elseif strcmp(Initialization, 'Squares')
%%% Squares initialization uses squares of sides half of the radius of
%%% preferred objects (no manual seeds are needed).

    if ~isempty(str2num(LayerNumString))
        if (str2num(LayerNumString)<1)
            LayerNumber = 4;
        else
            LayerNumber = str2num(LayerNumString);
        end
    end
    [InitialPhi] = createInitMLPhiSquares(ImageHeight, ImageWidth, PriorPhasefieldParameters, LayerNumber);
    ExtendedInitialPhi = zeros(ExtendedImageHeight, ExtendedImageWidth, LayerNumber) - 1;
    ExtendedInitialPhi(Maxd+1:Maxd+ImageHeight, Maxd+1:Maxd+ImageWidth, :) = (InitialPhi*2-1)*0.1;
    
else
%%% Neutral initialization sets the initial phase field to a Gaussian white
%%% noise with mean of alpha_f/lambda_f and a tiny variance (no manual
%%% seeds are needed).

    LayerNumber = str2num(LayerNumString);
    ExtendedInitialPhi = randn(ExtendedImageHeight,ExtendedImageWidth,LayerNumber) + ones(ExtendedImageHeight,ExtendedImageWidth,LayerNumber)*PriorPhasefieldParameters(1).alpha/PriorPhasefieldParameters(1).lambda;
end

ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);

%%% Check whether that figure is open. This checks all the figure handles
%%% for one whose handle is equal to the figure number for this module.
if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window.
    CPfigure(handles,'Image',ThisModuleFigureNumber);
    if handles.Current.SetBeingAnalyzed == handles.Current.StartingImageSet
        CPresizefigure(OrigImage,'TwoByTwo',ThisModuleFigureNumber)
    end
    %%% A subplot of the figure window is set to display the original image.
    h11 = subplot(2,2,1);
    CPimagesc(OrigImage,handles);
    title(['Input Image, cycle # ', num2str(handles.Current.SetBeingAnalyzed)]);
    
    %%% A subplot of the figure window is set to display the input objects
    %%% for correction
    h12 = subplot(2,2,2);
    im = CPlabel2rgb(handles,SegmentedGrayScaleObjects.*MaskSegmentedObject);
    CPimagesc(im,handles);
    title({'Masked Objects'; ['#layernum: ' num2str(LayerNumber)]});
end

drawnow;

%%% Geometric parameters of 'gas of circles' model to multi layered
%%% configuration
PriorPhasefieldParameters = repmat(PriorPhasefieldParameters, LayerNumber);

%%% Cut intensities over BG+4FG
%%% Using DataParameters fields
ExtendedImage(ExtendedImage>(DataParameters.muout+4*DataParameters.muin)) = DataParameters.muout+4*DataParameters.muin;

%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%

%%% Run the gradient descent algorithm of the multi layered gas of circles
%%% model
finalPhi = MLGOCSegmentationGM(PriorPhasefieldParameters, ExtendedImage, DataParameters, Kappa, ExtendedInitialPhi, OptimizationParameters);

%%% Set up the threshold to get the contour of foreground regions
Threshold = zeros(LayerNumber,1);

%%% Segmentation of image is get by thresholding
SegmentedLayers = zeros(size(finalPhi));
for ll=1:LayerNumber
    SegmentedLayers(:,:,ll) = finalPhi(:,:,ll) > Threshold(ll);
end

%%% Label multi dimensional binary image
SegmentedLayers = bwlabelml(SegmentedLayers);

%%%%%%%%%%%%%%%%%%%%%%
%%% POSTPROCESSING %%%
%%%%%%%%%%%%%%%%%%%%%%

%%% Remove embedded small objects
[EmbeddedRemovedObjects] = removeEmbeddedObjects( SegmentedLayers );

%%% Fill holes inside objects
for l=1:LayerNumber
    EmbeddedRemovedObjects(:,:,l) = imfill(EmbeddedRemovedObjects(:,:,l), 'holes');
end

%%% Remove objects area under pi*Radius*Radius/4 (the area of a circle with
%%% the radius half of the preferred
EmbeddedRemovedObjects = removeSmallObjects(EmbeddedRemovedObjects, Radius);

%%% Merge objects with Jaccard index over 0.8
MergedOverlappingObjects = mergeOverlappingObjects( EmbeddedRemovedObjects, 0.8 );

%%% From now we have the final objects satisfying size criteria and avoid
%%% degenerative cases; stored to find secondary objects
UneditedSegmentedImage = MergedOverlappingObjects;

%%% Remove overlapping objects
if OverlapJIThreshold>0.0
    [RemovedOverlappingObjects, RemovedOIdx] = removeOverlappingObjectsOverJIThreshold( MergedOverlappingObjects, OverlapJIThreshold );
else
    RemovedOverlappingObjects = MergedOverlappingObjects;
    RemovedOIdx = [];
end

%%% Remove overlapping objects
if strcmp(DiscardBorder, 'Yes')
    [RemovedOverlappingObjects, RemovedBIdx] = clearborderml(RemovedOverlappingObjects);
else
    RemovedBIdx = [];
end

FinalLabelMatrixImage = bwlabelml(RemovedOverlappingObjects);

%%% Calculating object perimeters
Perimeters = zeros(size(FinalLabelMatrixImage)) > 0;
for l = 1:LayerNumber
    Perimeters(:,:,l) = bwperim( FinalLabelMatrixImage(:,:,l)>0 );
end

FinalOutline = any(Perimeters, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SAVE DATA TO HANDLES STRUCTURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Saving Locations of the segmented objects.
handles.Measurements.(ObjectName).LocationFeatures = {'CenterX', 'CenterY'};
Centroids = zeros(max(FinalLabelMatrixImage(:)),2);

for l=1:LayerNumber
    tmp = FinalLabelMatrixImage(:,:,l);
    UniqueVales = unique(tmp);
    if ~isempty(UniqueVales) && length(UniqueVales)>1
        props = regionprops(tmp, 'Centroid');
        tmpCentroid = cat(1, props.Centroid);
        Centroids(UniqueVales(2:end)', :) = tmpCentroid(UniqueVales(2:end)', :);
    end
end
if isempty(Centroids)
    Centroids = [0 0];
end
handles.Measurements.(ObjectName).Location(handles.Current.SetBeingAnalyzed) = { Centroids };


%%% Saves the final phase field and thresholds to the handles structure.
handles.Pipeline.(['MLGOCFinalPhi' ObjectName]) = finalPhi;
handles.Pipeline.(['MLGOCThreshold' ObjectName]) = Threshold;

%%% Saves the segmented image, not edited for objects along the edges or
%%% for size, to the handles structure.
handles.Pipeline.(['UneditedSegmented' ObjectName]) = SegmentedLayers; % original segmentation without postprocessing

%%% Saves the segmented image, edited for small, embedded, to the
%%% handles structure.
handles.Pipeline.(['SmallRemovedSegmented' ObjectName]) = UneditedSegmentedImage; % contains overlapping objects

%%% Saves the final segmented label matrix image to the handles structure.
handles.Pipeline.(['Segmented' ObjectName]) = FinalLabelMatrixImage; % cleared from embedded and overlapping

%%% Saves the indices of removed objects; stored to discard secondary
%%% objects
handles.Pipeline.(['RemovedObjectIdx' ObjectName]) = [RemovedOIdx RemovedBIdx];

%%% Saves images to the handles structure so they can be saved to the hard
%%% drive, if the user requested.
if ~strcmpi(SaveOutlines,'Do not save')
    try
        handles.Pipeline.(SaveOutlines) = FinalOutline;
        handles.Pipeline.([SaveOutlines 'ML']) = Perimeters;
    catch e
        e.identifier
        error(['The object outlines were not calculated by the ', ModuleName, ' module, so these images were not saved to the handles structure. The Save Images module will therefore not function on these images. This is just for your information - image processing is still in progress, but the Save Images module will fail if you attempted to save these images.'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% DISPLAY RESULTS %%%
%%%%%%%%%%%%%%%%%%%%%%%
% drawnow

if any(findobj == ThisModuleFigureNumber)
    %%% Activates the appropriate figure window
    
    [ContourImage, Contours] = createContourImage(OrigImage, finalPhi, Threshold);
    
    %%% A subplot of the figure window is set to display the outlined image
    h21 = subplot(2,2,3);
    CPimagesc(ContourImage,handles);
    title('Segmentation');
    
    [ContourImageClear, ContoursClear] = createContourImage(OrigImage, (double(RemovedOverlappingObjects>0)*2-1)/2, Threshold);
    
    %%% A subplot of the figure window is set to display the sum of final
    %%% phase field layers
    h22 = subplot(2,2,4);
    CPimagesc(ContourImageClear, handles);
    title('Cleared');
end

end

%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

function extIm = extendImage(image, width, mode)
%EXTENDIM gives a padding to an image
%   EXTENDIM(IMAGE, WIDTH, MODE) returns matrix IMAGE padded with a border of
%   WIDTH. MODE can be 'mirror', 'torus' or a real number value.

[h, w, d] = size(image);
extIm = zeros(h+2*width, w+2*width, d);

switch mode
    case 'mirror' % the border of the image is mirrored
        for i=1:d
            
            extIm(1:width,1:width,i) = rot90(image(1:width,1:width,i),2); %nw
            extIm(1:width,width+w+1:end,i) =rot90(image(1:width,w-width+1:end,i),2); %ne
            extIm(width+h+1:end,1:width,i) = rot90(image(h-width+1:end,1:width,i),2); %dw
            extIm(width+h+1:end,width+w+1:end,i) = rot90(image(h-width+1:end,w-width+1:end,i),2); %de
            
            
            
            extIm(1:width,width+1:width+w,i) = flipud(image(1:width,:,i)); % top
            extIm(width+h+1:end,width+1:width+w,i) = flipud(image(h-width+1:end,:,i)); %bottom
            extIm(1+width:h+width,1:width,i) = fliplr(image(:,1:width,i)); % left
            extIm(1+width:h+width,width+w+1:end,i) = fliplr(image(:,w-width+1:end,i)); % right
        end
    case 'torus' % the border of the image is cut from the opposite border like a torus
        for i=1:d
            extIm(1:width,1:width,i) = image(h-width+1:end,w-width+1:end,i); %nw
            extIm(1:width,width+w+1:end,i) = image(h-width+1:end,1:width,i); %ne
            extIm(width+h+1:end,1:width,i) = image(1:width,w-width+1:end,i); %dw
            extIm(width+h+1:end,width+w+1:end,i) = image(1:width,1:width,i); %de
            
            extIm(1:width,width+1:width+w,i) = image(h-width+1:end,:,i); % top
            extIm(width+h+1:end,width+1:width+w,i) = image(1:width,:,i); %bottom
            extIm(1+width:h+width,1:width,i) = image(:,w-width+1:end,i); % left
            extIm(1+width:h+width,width+w+1:end,i) = image(:,1:width,i); % right
        end
    otherwise % padding with a constant value, usually by mean of the background of the image
        if isnumeric(mode)
            extIm = extIm + mode;
        end
end

% for i=1:d
%     extIm(width+1:width+h,width+1:width+w,i) = image(:,:,i); %original
% end

extIm(width+1:width+h,width+1:width+w,:) = image;

end

function initMLPhi = sortGrayscaleObjects2Layers(GrayscaleObjects, radius)
%%% SORTGRAYSCALEOBJECTS2LAYERS distributes the objects into saveral layers
%%% based on x and y coordinates

L = GrayscaleObjects;
num = max(GrayscaleObjects(:));

[h,w] = size(GrayscaleObjects);

layerNum = 1;

initMLPhi = zeros(h, w, layerNum);

SE = strel('disk',round(1.5*radius),0);

for i=1:num
    %     i
    objLayer = cutObject2Layer(L,i);
    
    l=1; 
    minOverlap = h*w;
    %     minIndex = 1;
    ready = 0;
    while l<=layerNum && ~ready
        overlap = numel( find( initMLPhi(:,:,l) & imdilate(objLayer,SE) ) );
        if overlap == 0
            %             'no overlap';
            initMLPhi(:,:,l) = initMLPhi(:,:,l) | objLayer;
            ready = 1;
        else
            if overlap<minOverlap
                %                 'smaller overlap';
                minOverlap = overlap;
                %                 minIndex = l;
            end
        end
        
        l = l+1;
    end
    
    if ready == 0
        layerNum = layerNum+1;
        initMLPhi(:,:,layerNum) = objLayer;
    end
    
end

end

function initMLPhi = sortGrayscaleObjects2LayersColoring(GrayscaleObjects, ncolor)
%%% SORTGRAYSCALEOBJECTS2LAYERSCOLORING distributes the objects into saveral
%%% layers based on map coloring algorithm

L = GrayscaleObjects;

[h,w] = size(GrayscaleObjects);

layerNum = ncolor;

initMLPhi = zeros(h, w, layerNum);

[colors, ~] = patchColorSelect(L, ncolor);

for ll = 1:layerNum
    objectIndicesOfLayer = find(colors==ll);
    tempLayer = zeros(h, w);
    for oi = objectIndicesOfLayer
        tempLayer(L==oi) = 1;
    end
    initMLPhi(:,:,ll) = tempLayer;
end

end

function [colors, pc] = patchColorSelect( GrayscaleObjects, ncolor )

props = regionprops(GrayscaleObjects, 'Centroid');
centers = cat(1, props.Centroid);
iNotEmpty = zeros(length(props),1) +1;


for objectNum = length(props):-1:1
    centers(objectNum, :);
    if any(isnan(centers(objectNum,:)))
        centers(objectNum,:) = [];
        iNotEmpty(objectNum) = 0;
    end
end

if nargin < 2
    ncolor = 7;
end
ncolor = min(ncolor,7);
[~,ic] = sort(sum(abs(centers).^2,2));
n = size(centers,1);
colors(ic) = mod(0:n-1,ncolor)+1;
pc = {'b.','g.','r.','c.','m.','y.','k.'};
pc = pc(1:ncolor);
for k1=ncolor+1:n
    k0 = ic(1:k1-1);
    k = ic(k1);
    cd = sum(abs(centers(k0,:) - centers(k+zeros(k1-1,1),:)).^2,2);
    s = 0;
    cc = colors(k0);
    for c = 1:ncolor
        id = min(cd(cc==c));
        if id > s
            s = id;
            is = c;
        end
    end
    colors(k) = is;
    kc = find(cd'==s & cc==is, 1);
    kk = k0(kc);
    k0(kc) = k;
    k = kk;
    kd = sum(abs(centers(k0,:) - centers(k+zeros(numel(k0),1),:)).^2,2);
    ks = 0;
    cc = colors(k0);
    for c = 1:ncolor
        id = min(kd(cc==c));
        if id > ks
            ks = id;
            ik = c;
        end
    end
    if ks > s
        colors(k) = ik;
    end
end
end

function [L, num] = splitObjects(L, num, radius)
%%% SPLITOBJECTS slices the objects over a specific size

%     'splitobjects started..'
maxNum = num;
radius = uint32(radius);
for i=1:num
    objLayer = uint32(cutObject2Layer(L,i));
    
    props = regionprops(objLayer, 'BoundingBox');
    
    wo = uint32(props.BoundingBox(3));
    ho = uint32(props.BoundingBox(4));
    
    x0 = uint32(props.BoundingBox(1));
    y0 = uint32(props.BoundingBox(2));
    
    numh = uint32(floor(ho/(2*radius*1.0)));
    numw = uint32(floor(wo/(2*radius*1.0)));
    
    if numh>numw
        extraLabelNum = numh-1;
        if extraLabelNum>0
            maxNum = maxNum + extraLabelNum;
            for s = 2:numh-1
                objLayer(y0+(s-1)*2*(radius):y0+(s)*2*(radius)-1,x0:x0+wo) = objLayer(y0+(s-1)*2*(radius):y0+(s)*2*(radius)-1,x0:x0+wo)*(maxNum - (extraLabelNum+s-1));
            end
            objLayer(y0+(numh-1)*2*radius:y0+ho,x0:x0+wo) = objLayer(y0+(numh-1)*2*radius:y0+ho,x0:x0+wo)*(maxNum);
        end
        
        
    else
        extraLabelNum = numw-1;
        if extraLabelNum>0
            maxNum = maxNum + extraLabelNum;
            
            for s = 2:numh-1
                objLayer(y0:y0+ho,x0+(s-1)*2*radius:x0+(s)*2*radius-1) = objLayer(y0:y0+ho,x0+(s-1)*2*radius:x0+(s)*2*radius-1)*(maxNum - (extraLabelNum+s-1));
            end
            objLayer(y0:y0+ho,x0+(numw-1)*2*radius:x0+wo) = objLayer(y0:y0+ho, x0+(numw-1)*2*radius:x0+wo)*(maxNum);
        end
    end
   
    for n = maxNum-extraLabelNum+1:maxNum
        lIndices = objLayer==n;
        L(lIndices) = n;
    end
    
end
num = maxNum;

end

function [objLayer] = cutObject2Layer( image, n )
%%% CUTOBJECT2LAYER gets an object with a specific label from a labeled
%%% image

[H, W] = size(image);
[X, Y] = find( image==n );

objLayer = uint8(zeros( H, W ));

N = size( X, 1 );
for i=1:N
    objLayer( X(i), Y(i) ) = 1;
end
end

function [contourImage, contours] = createContourImage(inputImage, MLPhi, threshold)
%%% CREATECONTOURIMAGE generates a colored contour image from a phase field
%%% segmentation

colors = [	1 0 0;
    0 1 0;
    0 0 1;
    1 1 0;
    1 0 1;
    0 1 1;
    0.5 0 0;
    0 0.5 0;
    0 0 0.5;
    0.5 0.5 0;
    0.5 0 0.5;
    0 0.5 0.5;
    0.5 0.5 0.5];

[h, w, layerNum] = size(MLPhi);

imInfo = whos('inputImage');

[hi, wi] = size(inputImage);

if h~=hi || w~=wi
    error('Image and phi must have same size!')
end

%     SE = strel('disk',1);

resPhi = zeros(size(MLPhi));
edgesPhi = zeros(size(MLPhi))>0;

%     figure;
for i=1:layerNum
    resPhi(:,:,i) = double(MLPhi(:,:,i) > threshold(i));
    
    %       edgesPhi(:,:,i) = xor(imerode(resPhi(:,:,i),SE),resPhi(:,:,i));
    edgesPhi(:,:,i) = bwperim(resPhi(:,:,i));
end

%     pause;

sumLayers = sum(edgesPhi, 3);
sumLayersIndex = sumLayers>0;
if strcmp(imInfo.class, 'uint8')
    contourLayer = double(inputImage)/2^8;
elseif strcmp(imInfo.class, 'uint16')
    contourLayer = double(inputImage)/2^16;
else
    contourLayer = inputImage;
end

% histogram stretch
contourLayer = imadjust(contourLayer, [min(contourLayer(:)) max(contourLayer(:))], [0.1 0.9]);

contourLayer(sumLayersIndex) = 0.0;

contourImage = zeros(h,w,3);

contourImage(:,:,1) = contourLayer;
contourImage(:,:,2) = contourLayer;
contourImage(:,:,3) = contourLayer;

sumLayers(sumLayers==0) = 1;

for i = 1:layerNum
    contourImage(:,:,1) = contourImage(:,:,1) + colors(mod(i, size(colors,1)),1)*(edgesPhi(:,:,i)./sumLayers);
    contourImage(:,:,2) = contourImage(:,:,2) + colors(mod(i, size(colors,1)),2)*(edgesPhi(:,:,i)./sumLayers);
    contourImage(:,:,3) = contourImage(:,:,3) + colors(mod(i, size(colors,1)),3)*(edgesPhi(:,:,i)./sumLayers);
end


contours = edgesPhi>0;
%     contourImage(isnan(contourImage)) = contourLayer;

end

function [initMLPhi] = createInitMLPhiSquares(Nx, Ny, parameters, layerNum)
%%% CREATEINITMLPHISQUARES creates initial phase field by covering the
%%% whole image space by squares (no initial objects are needed)

d = parameters(1).d;
d = int32(round(d/2));

if layerNum < 2
    layerNum = 2;
end

% bgSquare = zeros(2*d, 2*d) - 1;
bgSquare = zeros(2*d, 2*d) ;
fgSquare = zeros(2*d, 2*d) + 1;

if layerNum==2
    
    initMLPhi = zeros(Nx, Ny, layerNum);
    
    tile = [bgSquare fgSquare; fgSquare bgSquare];
    [hx,wx] = size(tile);
    
    xSize = round( Nx / hx );
    ySize = round( Ny / wx );
    
    layer1 = repmat(tile, xSize+1, ySize+1);
    
    initMLPhi(:,:,1) = layer1(1:Nx, 1:Ny);
    
    tile = [fgSquare bgSquare; bgSquare fgSquare];
    
    [hx,wx] = size(tile);
    
    xSize = round( Nx / hx );
    ySize = round( Ny / wx );
    
    layer2 = repmat(tile, xSize+1, ySize+1);
    
    initMLPhi(:,:,2) = layer2(1:Nx, 1:Ny);
    
else
    
    %     tile = zeros(2*d, 2*d*layerNum) - 1;
    tile = zeros(2*d, 2*d*layerNum) ;
    tile(:, 1:2*d) = 1;
    
    xSize = round( Nx/(2*d*layerNum) );
    
    bigTileRow1 = repmat( tile, 1, xSize+2);
    
    %     bigTileRow2 = zeros(size(bigTileRow1)) - 1;
    bigTileRow2 = zeros(size(bigTileRow1)) ;
    
    bigTileRow2(:, 3*d+1:end) = bigTileRow1( :, 1:end-3*d );
    
    ySize = round( Ny/(4*d) );
    
    patternAll = repmat( [bigTileRow1; bigTileRow2], ySize+1, 1);
    
    initMLPhi = zeros(Nx, Ny, layerNum);
    
    for i=1:layerNum
        %         Nx+(i-1)*2*d
        initMLPhi(:,:,i) = patternAll(1:Nx, 1+(i-1)*2*d:Ny+(i-1)*2*d);
    end
    
end

end

function InitialCircularPhi = createInitMLPhiCircles(InitialPhi, Radius)
%%% CREATEINITMLPHICIRCLES creates initial phase field by replacing each
%%% object with a circle at the centroid of the original object

InitialCircularPhi = zeros(size(InitialPhi));

[x, y] = meshgrid(1:2*Radius, 1:2*Radius);
    
circle = uint8(sqrt((x-Radius-0.5).^2+(y-Radius-0.5).^2) <= Radius);

LabeledInitialPhi = bwlabelml(InitialPhi);

props = regionprops(LabeledInitialPhi, 'Centroid');
% Centroids = cat(1, props.Centroid);

for i=1:length(props)
    
%     [x0, y0, l0] = props(i).Centroid;
    x0 = round(props(i).Centroid(1));
    y0 = round(props(i).Centroid(2));
    l0 = round(props(i).Centroid(3));
    InitialCircularPhi(x0-Radius+1:x0+Radius, y0-Radius+1:y0+Radius, l0) = circle*InitialPhi(x0, y0, l0);
    
end

end

function [OverlappingObjectIndexes, OOL] = getOverlappingObjectIndexes( SegmentedLayers )
%%% GETOVERLAPPINGOBJECTINDEXES collects the indices of objects that has
%%% any overlap with each other

%%% Storing where are overlapping objects
% OverlapMap = zeros(ImageHeight, ImageWidth);
% BgIdx = find(sum(SegmentedLayers>0, 3) == 0);
% SingleIdx = find(sum(SegmentedLayers>0, 3) == 1);
% OverlapIdx = find(sum(SegmentedLayers>0, 3) > 1);

[ImageHeight, ImageWidth, LayerNumber] = size(SegmentedLayers);

OverlapMap = zeros(ImageHeight, ImageWidth);
BgIdx = (sum(SegmentedLayers>0, 3) == 0);
SingleIdx = (sum(SegmentedLayers>0, 3) == 1);
OverlapIdx = (sum(SegmentedLayers>0, 3) > 1);

OverlapMap(BgIdx) = 0;
OverlapMap(SingleIdx) = 1;
OverlapMap(OverlapIdx) = 2;

OverlappingObjectIndexes = [];
OOL = [];
indDOOI = 0;

SE = strel('disk',1,4);

for ll=1:LayerNumber
    
    UniqueValues = unique(SegmentedLayers(:,:,ll));
    
    if ~isempty(UniqueValues) && length(UniqueValues)>1
        for Object = UniqueValues(2:end)'
            
            hitTest = (imerode(SegmentedLayers(:,:,ll)==Object, SE)) & (OverlapMap==2);
            
            %%% check the object has overlap area
            if sum(hitTest(:))>0
                
                indDOOI = indDOOI+1;
                OverlappingObjectIndexes(indDOOI) = Object;
                OOL(indDOOI) = ll;
                
            end
        end
    end
    
end

end

function RemovedEmbeddedObjecs = removeEmbeddedObjects( SegmentedLayers )
%%% REMOVEEMBEDDEDOBJECTS removes each 'small' objects that are surrounded
%%% by a bigger object

display('Removing embedded objects...');

LayerNumber = size(SegmentedLayers, 3);

RemovedEmbeddedObjecs = SegmentedLayers;
DiscardEmbeddedObjIdx = [];
indDEOI = 0;

% if strcmp(DiscardEmbedded, 'Yes')

for ll=1:LayerNumber
    
    UniqueValues = unique(RemovedEmbeddedObjecs(:,:,ll));
    
    if ~isempty(UniqueValues) && length(UniqueValues)>1
        
        for Object = UniqueValues(2:end)'
            
            for lll = 1:LayerNumber
                if lll==ll
                    continue;
                else
                    if (or( RemovedEmbeddedObjecs(:,:,ll)==Object, RemovedEmbeddedObjecs(:,:,lll)>0 )) == (RemovedEmbeddedObjecs(:,:,lll)>0)
                        indDEOI = indDEOI + 1;
                        DiscardEmbeddedObjIdx(indDEOI) = Object;
                    end
                end
            end
        end
    end
    
end

for indO=1:length(DiscardEmbeddedObjIdx)
    RemovedEmbeddedObjecs(RemovedEmbeddedObjecs==DiscardEmbeddedObjIdx(indO)) = 0;
end

% for ll=1:LayerNumber
%     
%     SegmentedLayer = bwlabel(RemovedEmbeddedObjecs(:,:,ll));
%     if ll==1
%         %         MaxId = max( reshape(SegmentedLayer,1,[]) );
%         RemovedEmbeddedObjecs(:,:,ll) = SegmentedLayer;
%     else
%         MaxId = max( reshape(RemovedEmbeddedObjecs(:,:,ll-1), 1, []) );
%         SegmentedLayer(SegmentedLayer>0) = SegmentedLayer(SegmentedLayer>0) + MaxId;
%         RemovedEmbeddedObjecs(:,:,ll) = SegmentedLayer;
%     end
% end

[RemovedEmbeddedObjecs, num] = bwlabelml(RemovedEmbeddedObjecs);
% end

if isempty(RemovedEmbeddedObjecs)
    RemovedEmbeddedObjecs = zeros(ImageHeight, ImageWidth);
end

display([num2str(length(DiscardEmbeddedObjIdx)) ' embedded objects were removed.']);

end

function MergedOverlappingObjects = mergeOverlappingObjects( SegmentedLayers, JIThreshold)
%%% MERGEOVERLAPPINGOBJECTS checks all overlapping object pairs and merges
%%% them if the Jaccard index of them is higher then a threshold

display('Merging fully overlapping objects...');

[OverlappingObjectIndexes, OOL] = getOverlappingObjectIndexes( SegmentedLayers );

RemovedIdxList = [];

for ooiInd = 1:length(OverlappingObjectIndexes)
    ooi = OverlappingObjectIndexes(ooiInd);
    
    if ~ismember(ooi, RemovedIdxList)
        % to avoid double checking we start second indexing from ooiInd
        for oojInd = ooiInd:length(OverlappingObjectIndexes)
            ooj = OverlappingObjectIndexes(oojInd);
            
            if ~ismember(oojInd, RemovedIdxList) && ooi~=ooj
                ObjectI = SegmentedLayers(:,:,OOL(ooiInd)) == ooi;
                ObjectJ = SegmentedLayers(:,:,OOL(oojInd)) == ooj;
                
                if jaccardCoefficient(ObjectI, ObjectJ)>JIThreshold
                    MergedObject = any( cat(3, ObjectI, ObjectJ) ,3);
                    SegmentedLayers(:,:,OOL(ooiInd)) = SegmentedLayers(:,:,OOL(ooiInd)) - ObjectI*ooi + MergedObject*ooi;
                    SegmentedLayers(:,:,OOL(oojInd)) = SegmentedLayers(:,:,OOL(oojInd)) - ObjectJ*ooj;
                    RemovedIdxList = [RemovedIdxList ooj];
                end
            end
            
        end
    end
end

[MergedOverlappingObjects, num] = bwlabelml(SegmentedLayers);

display([num2str(length(RemovedIdxList)) ' object(s) were merged (and discarded).']);

end

function [RemovedOverlappingObjects, RemovedIdxList] = removeOverlappingObjectsOverJIThreshold( SegmentedLayers, OverlapJIThreshold )
%%% REMOVEOVERLAPPINGOBJECTSOVERJITHRESHOLD checks every pair of
%%% overlapping objects, and discards the bigger object, if the Jaccard
%%% index of the overlapping area and the smaller object is over
%%% OverlapJIThreshold

fprintf('Removing objects with Jaccard index over %0.2f\n', OverlapJIThreshold);

[OverlappingObjectIndexes, OOL] = getOverlappingObjectIndexes( SegmentedLayers );

RemovedIdxList = [];

for ooiInd = 1:length(OverlappingObjectIndexes)
    ooi = OverlappingObjectIndexes(ooiInd);
    ObjectI = SegmentedLayers(:,:,OOL(ooiInd)) == ooi;
    PropsI = regionprops(ObjectI, 'Area', 'PixelIdxList','Centroid');
    
    if ~ismember(ooi, RemovedIdxList)
        % to avoid double checking we start second indexing from ooiInd
        for oojInd = ooiInd:length(OverlappingObjectIndexes)
            ooj = OverlappingObjectIndexes(oojInd);
            
            if ~ismember(oojInd, RemovedIdxList) && ooi~=ooj
                
                ObjectJ = SegmentedLayers(:,:,OOL(oojInd)) == ooj;
                PropsJ = regionprops(ObjectJ, 'Area', 'PixelIdxList','Centroid');
                
                if (PropsI.Area>PropsJ.Area)
%                     BiggerObject = ObjectI;
                    PropsBiggerObject = PropsI;
                    IdxBiggerObject = ooi;
%                     SmallerObject = ObjectJ;
                    PropsSmallerObject = PropsJ;
%                     IdxSmallerObject = ooj;
                else
%                     BiggerObject = ObjectJ;
                    PropsBiggerObject = PropsJ;
                    IdxBiggerObject = ooj;
%                     SmallerObject = ObjectI;
                    PropsSmallerObject = PropsI;
%                     IdxSmallerObject = ooi;
                end
                
                if ismember(IdxBiggerObject, RemovedIdxList)
                    continue;
                end
                
                JI = jaccardCoefficientByPixelIdxList(PropsSmallerObject.PixelIdxList, intersect(PropsSmallerObject.PixelIdxList,PropsBiggerObject.PixelIdxList));

                if JI>OverlapJIThreshold % && 2*numel(PropsSmallerObject.PixelIdxList)<numel(PropsBiggerObject.PixelIdxList)
%                     MergedObject = any( cat(3, ObjectI, ObjectJ) ,3);
%                     SegmentedLayers(:,:,OOL(ooiInd)) = SegmentedLayers(:,:,OOL(ooiInd)) - ObjectI*ooi + MergedObject*ooi;
%                     SegmentedLayers(:,:,OOL(oojInd)) = SegmentedLayers(:,:,OOL(oojInd)) - ObjectJ*ooj;
                    RemovedIdxList = [RemovedIdxList IdxBiggerObject];
                    
%                     figure(11);
%                     imagesc( (BiggerObject>0)*-1 + (SmallerObject>0)*1 );
%                     
%                     text(PropsBiggerObject.Centroid(1), PropsBiggerObject.Centroid(2), sprintf('BID: %d, area: %d',IdxBiggerObject,PropsBiggerObject.Area),'EdgeColor',[0 0 1]);
%                     text(PropsSmallerObject.Centroid(1), PropsSmallerObject.Centroid(2), sprintf('BID: %d, area: %d',IdxSmallerObject,PropsSmallerObject.Area),'EdgeColor',[1 0 0]);
%                     
%                     title(sprintf('%d vs %d: JI=%0.2f',IdxBiggerObject,IdxSmallerObject,JI));
                    
                end
                
                
            end
            
        end
    end
end

for indO = 1:length(RemovedIdxList)
    SegmentedLayers( SegmentedLayers == RemovedIdxList(indO) ) = 0;
end

if isempty(SegmentedLayers)
    SegmentedLayers = uint8(zeros(ImageHeight, ImageWidth));
end

[RemovedOverlappingObjects, num] = bwlabelml(SegmentedLayers);

display([num2str(length(RemovedIdxList)) ' object(s) were discarded.']);

end

function [RemovedOverlappingObjects, RemovedIdx] = removeOverlappingObjects( SegmentedLayers )
%%% REMOVEOVERLAPPINGOBJECTS discards the objects that have overlapping to
%%% each other and stores their indices

display('Removing overlapping objects...');

[ImageHeight, ImageWidth, LayerNumber] = size( SegmentedLayers );

% SE = strel('disk', 1, 4);

[OverlappingObjectIndexes, OOL] = getOverlappingObjectIndexes( SegmentedLayers );

for indO = 1:length(OverlappingObjectIndexes)
    SegmentedLayers( SegmentedLayers == OverlappingObjectIndexes(indO) ) = 0;
end

if isempty(SegmentedLayers)
    SegmentedLayers = uint8(zeros(ImageHeight, ImageWidth));
end

RemovedOverlappingObjects = SegmentedLayers;

RemovedIdx = OverlappingObjectIndexes;

display([num2str(length(RemovedIdx)) ' overlapping objects were removed.']);

end

function RemovedSmallObjects = removeSmallObjects(RemovedOverlappingObjects, MinSize)
%%% REMOVESMALLOBJECTS discards all objects below a given size

display('Removing small objects...');

[ImageHeight, ImageWidth, LayerNumber] = size(RemovedOverlappingObjects);
RemovedSmallObjects = zeros(ImageHeight, ImageWidth, LayerNumber);

for l = 1:LayerNumber
    tmp = RemovedOverlappingObjects(:,:,l);
    %%% Get diameters of objects and calculate the interval
    %%% that contains 90% of the objects
    props = regionprops(RemovedOverlappingObjects(:,:,l),'EquivDiameter');
    Diameters = [0;cat(1,props.EquivDiameter)];
%     SortedDiameters = sort(Diameters);
    %     NbrInTails = max(round(0.05*length(Diameters)),1);
    %     Lower90Limit = SortedDiameters(NbrInTails);
    %     Upper90Limit = SortedDiameters(end-NbrInTails+1);
    
    %%% Locate objects with diameter outside the specified range
    
    %%% Create image with object intensity equal to the diameter
    DiameterMap = Diameters(tmp+1);
    %%% Remove objects that are too small
    tmp(DiameterMap < MinSize) = 0;
    
    RemovedSmallObjects(:,:,l) = tmp;
end

[RemovedSmallObjects, num] = bwlabelml(RemovedSmallObjects);

display([num2str(numel(unique(RemovedOverlappingObjects(:)))-numel(unique(RemovedSmallObjects(:)))) ' small objects were removed.']);

end

function [RemovedBorderObjects, RemovedIdx] = clearborderml(SegmentedLayers)
%%% CLEARBORDERML discards the objects touching the image border and stores
%%% their indices

display('Removing border objects...');

[ImageHeight, ImageWidth, LayerNumber] = size(SegmentedLayers);
RemovedBorderObjects = zeros(ImageHeight, ImageWidth, LayerNumber);

for ll=1:LayerNumber
    ActLayer = SegmentedLayers(:,:,ll);
    map = 0:max(ActLayer(:));
    map(ActLayer(1,:)+1) = 0;
    map(ActLayer(end,:)+1) = 0;
    map(ActLayer(:,1)+1) = 0;
    map(ActLayer(:,end)+1) = 0;
    ActLayer = map(ActLayer+1);
    RemovedBorderObjects(:,:,ll) = ActLayer;
end

% [RemovedBorderObjects, num] = bwlabelml(RemovedBorderObjects);

RemovedIdx = unique(SegmentedLayers-RemovedBorderObjects)';
RemovedIdx = RemovedIdx(RemovedIdx~=0);

display([num2str(numel(RemovedIdx)) ' border objects were removed.']);

end

function [jaccardIdx,jaccardDist] = jaccardCoefficient(img_Orig,img_Seg)
% Jaccard index and distance co-efficient of segmented and ground truth
% image
% Usage: [index,distance(JC)] = jaccard_coefficient(Orig_Image,Seg_Image);

% Check for logical image
if ~islogical(img_Orig)
    error('Image must be in logical format');
end
if ~islogical(img_Seg)
    error('Image must be in logical format');
end

% Find the intersection of the two images
inter_image = img_Orig & img_Seg;

% Find the union of the two images
union_image = img_Orig | img_Seg;

if sum(union_image(:))~=0
    jaccardIdx = sum(inter_image(:))/sum(union_image(:));
else
    jaccardIdx = 0;
end
% Jaccard distance = 1 - jaccardindex;
jaccardDist = 1 - jaccardIdx;

end

function [jaccardIdx,jaccardDist] = jaccardCoefficientByPixelIdxList(PixelIdx1, PixelIdx2)

jaccardIdx = numel( intersect(PixelIdx1,PixelIdx2) ) / numel( union(PixelIdx1, PixelIdx2) );
jaccardDist = 1 - jaccardIdx;

end

function [LML, num] = bwlabelml(L)
% BWLABELML is an extension of BWLABEL function to multi dimensional binary
% images.

[ImageHeight, ImageWidth, LayerNumber] = size(L);
LML = zeros(ImageHeight, ImageWidth, LayerNumber);
MaxId = 0;

for ll=1:LayerNumber
    
    SegmentedLayer = bwlabel( L(:,:,ll) > 0 );
    if ll==1
        LML(:,:,ll) = SegmentedLayer;
    else
        SegmentedLayer(SegmentedLayer>0) = SegmentedLayer(SegmentedLayer>0) + MaxId;
        LML(:,:,ll) = SegmentedLayer;
    end
    MaxId = max(LML(:));
    
end

num = MaxId;

end