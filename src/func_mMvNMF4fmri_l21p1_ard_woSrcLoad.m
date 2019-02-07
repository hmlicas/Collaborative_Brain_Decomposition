function [U, V] = func_mMvNMF4fmri_l21p1_ard_woSrcLoad(sbjData,prepData,outDir,resId,initName,sbjNum,K,alphaS21,alphaL,spaR,vxI,ard,eta,iterNum,calcGrp,parforOn)

%disp(['# of input: ',num2str(nargin)]);
if nargin~=16
    error('16 input parameters are needed');
end

% if isdeployed
%     sbjNum = str2double(sbjNum);
%     K = str2double(K);
%     alphaS21 = str2double(alphaS21);
%     alphaL = str2double(alphaL);
%     spaR = str2double(spaR);
% 	vxI = str2double(vxI);
% 	ard = str2double(ard);
% 	eta = str2double(eta);
% 	iterNum = str2double(iterNum);
%     calcGrp = str2double(calcGrp);
%     parforOn = str2double(parforOn);
% end

alphaS1 = alphaS21;
alphaS21 = round(alphaS21*sbjNum);

% parameter setting
options = [];
options.maxIter = 100;
options.error = 1e-4;
options.nRepeat = 1;
options.minIter = 30;
options.meanFitRatio = 0.1;
options.rounds = iterNum;
options.NormW = 1;
options.eta = eta;

options.alphaS21 = alphaS21;
options.alphaL = alphaL;
options.vxlInfo = vxI;
options.spaR = spaR;
if ard>0
	options.ardUsed = ard;
end

% output name
resDir = [outDir,filesep,resId,'_sbj',num2str(sbjNum),'_comp',num2str(K),...
		  '_alphaS21_',num2str(alphaS21),'_alphaL',num2str(alphaL),...
		  '_vxInfo',num2str(vxI),'_ard',num2str(ard),'_eta',num2str(eta)];
if ~exist(resDir,'dir')
    mkdir(resDir);
end

% load data
numUsed = sbjNum;
numUsed = min(numUsed,length(sbjData));
sbjData = sbjData(1:numUsed);

% load preparation, containing gNb
load(prepData);

% load initialization: initV
load(initName,'initV');

%
vNum = length(gNb);
W = cell(numUsed,1);
Wt = cell(numUsed,1);
D = cell(size(W));
L = cell(size(W));
Dt = cell(size(Wt));
Lt = cell(size(Wt));
U = cell(numUsed,1);
V = cell(numUsed,1);

disp('preprocess...');
for si=1:numUsed
    fprintf('  sbj%d->',si);
    if mod(si,15)==0
        fprintf('\n');
    end

    origSbjData = sbjData{si};
    origSbjData = dataPrepro(origSbjData,'vp','vmax');
    sbjData{si} = origSbjData;

    nanSbj = isnan(origSbjData);
    if sum(nanSbj(:))>0
        disp([' nan exists: ','sbj',num2str(si)]);
        return;
    end

    % construct the spatial affinity graph
    if vxI==0
        if si==1
        	tmpW = sparse(vNum,vNum);
            for vi=1:vNum
                for ni=1:length(gNb{vi})
			        nei = gNb{vi}(ni);
			        tmpW(vi,nei) = 1;
			        tmpW(nei,vi) = 1;
                end
            end
            W{si} = tmpW;
        else
            W{si} = W{1};
        end
    else
    	tmpW = sparse(vNum,vNum);
        for vi=1:vNum
            for ni=1:length(gNb{vi})
                nei = gNb{vi}(ni);
                if vi<nei
                    corrVal = (1+corr(origSbjData(:,vi),origSbjData(:,nei)))/2;
                    if isnan(corrVal)
			            corrVal = 0;
                    end
			        tmpW(vi,nei) = corrVal;
			        tmpW(nei,vi) = corrVal;
	            else
	                continue;
                end
            end
        end
		W{si} = tmpW;
    end

    % temporal affinity matrix
	if isfield(options,'alphaLT')
		t = size(origSbjData,1);
		Wt{si} = zeros(t,t);
		if isfield(options,'timR')
			tNei = timR;
		else
			tNei = 1;
		end
		for tni=1:tNei
			Wt{si} = Wt{si} + diag(ones(t-tni,1),tni) + diag(ones(t-tni,1),-tni);
		end
	else
		Wt{si} = [];
	end

	[mFea,nSmp] = size(sbjData{si});
    
    if isfield(options,'alphaL')
		DCol = full(sum(W{si},2));
		D{si} = spdiags(DCol,0,nSmp,nSmp);
		L{si} = D{si} - W{si};
		if isfield(options,'NormW') && options.NormW
			D_mhalf = spdiags(DCol.^-0.5,0,nSmp,nSmp);
			L{si} = D_mhalf*L{si}*D_mhalf * options.alphaL;
			W{si} = D_mhalf*W{si}*D_mhalf * options.alphaL;
			D{si} = D_mhalf*D{si}*D_mhalf * options.alphaL;
		end
	else
		D{si} = [];
		L{si} = [];
    end
    
    if isfield(options,'alphaLT')
		DCol = full(sum(Wt{si},2));
		Dt{si} = spdiags(DCol,0,mFea,mFea);
		Lt{si} = Dt{si} - Wt{si};
		if isfield(options,'NormW') && options.NormW
			D_mhalf = spdiags(DCol.^-0.5,0,mFea,mFea);
			Lt{si} = D_mhalf*Lt{si}*D_mhalf * options.alphaLT;
			Wt{si} = D_mhalf*Wt{si}*D_mhalf * options.alphaLT;
			Dt{si} = D_mhalf*Dt{si}*D_mhalf * options.alphaLT;
		end
	else
		Dt{si} = [];
		Lt{si} = [];
    end
    
	% initialization (old)
	V{si} = initV;
	miv = max(V{si});
	trimInd = V{si}./max(repmat(miv,size(V{si},1),1),eps) < 1e-2;
	V{si}(trimInd) = 0;
    
    U_ = [];
	U{si} = backNMF_u(sbjData{si}, K, options, U_, V{si});

	% initialization (updated)
%     U_ = [];
% 	dualOpt = [];
%     dualOpt.maxIter = 500;
%     dualOpt.error = 1e-4;
%     dualOpt.nRepeat = 1;
%     dualOpt.minIter = 100;
%     dualOpt.meanFitRatio = 0.1;
%     dualOpt.NormW = 1;
%     dualOpt.S1 = alphaS1;
%     dualOpt.L = alphaL;
%     [U{si},~] = backNMF_u(sbjData{si},K,dualOpt,U_,V{si});
%     
%     [~,V{si}] = mNMF_sp_v(sbjData{si},K,W{si},dualOpt,U{si},V{si});
%     miv = max(V{si});
%     trimInd = V{si}./max(repmat(miv,size(V{si},1),1),eps) < 1e-2;
%     V{si}(trimInd) = 0;
end
fprintf('\n');
initUvName = [resDir,filesep,'init_UV.mat'];
save(initUvName,'U','V','-v7.3');

% decomposition l21
if parforOn==0
    [U, V] = mMultiNMF_l21p1_ard(sbjData, W, D, L, Wt, Dt, Lt, U, V, options); 
elseif parforOn==1
    [U, V] = mMultiNMF_l21p1_ard_parfor(sbjData, W, D, L, Wt, Dt, Lt, U, V, options); 
end

% finalize the results
finalUvName = [resDir,filesep,'final_UV.mat'];
save(finalUvName,'U','V','-v7.3');

% compute the V_centroid
V_centroid = zeros(size(V{1}));
for vi=1:length(V)
    V_centroid = V_centroid + V{vi};
end
V_centroid = V_centroid / length(V);

outName_cen = [resDir,filesep,'res_cen.mat'];
save(outName_cen,'V_centroid','-v7.3');

if calcGrp==1
    disp('calculate time course based on group spatial maps (for comparison use later)...');
    
    grpOpt = [];
    grpOpt.maxIter = 500;
    grpOpt.error = 1e-4;
    grpOpt.nRepeat = 1;
    grpOpt.minIter = 100;
    grpOpt.meanFitRatio = 0.1;
    grpOpt.NormW = 1;
    
    bSbjNum = length(V);
    U = cell(bSbjNum,1);
    V = cell(bSbjNum,1);
    for si=1:bSbjNum
        [U{si},~] = backNMF_u(sbjData{si}, K, grpOpt, [], V_centroid);
        V{si} = V_centroid;
    end

    outName_grp = [resDir,filesep,'grp_UV.mat'];
    save(outName_grp,'U','V','-v7.3');
end

% if isdeployed
%     exit;
% else
%     disp('Done!');
% end


