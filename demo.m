% Demo function for DTVNL1 Kornecker Compressive Sensing recovery in ICIP paper 
% 
% Thuong Nguyen Canh, Khanh Quoc Dinh, and Byeungwoo Jeon, "Detail-preserving compressive
% sensing recovery based on cartoon texture image decomposition", IEEE Inter. Conf. Image.
% Process. (ICIP), Paris, 2014. 
% 
% Please CITE this paper if you are using this source code. 
%
% DTVNL1: Decomposition based Total Variation with Spatial Nonlocal Regularization 
% Copyright (C) 2015  Thuong Nguyen Canh
% https://sites.google.com/site/ngcthuong 
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


close all;  clear all; % clc;
path(path,genpath(pwd));

%% General setting
par.imgSize     = 256;
par.sparsity    = [0.2];
testIm          = [1 ];
par.testIm      = testIm; 
par.recMode     = 'DTVNL1'; % 
par.filterMode  = 'BM3D';   % 'No', 'NLM', 'BM3D'
par.sigma       = 0.02;    
par.showPSNR    = true;  % true/false
par.nbrLoop     = 1;
par.maxIter     = 25; 
par.qid_index   = [1 2 3  5:19];
par.testIm      = 5;  
custorm_note    = ['_20141127_iterMax' num2str(par.maxIter)]; 

for imgId = 1:1:length(testIm)
    [imgOrg, imgName] = testImage(par.imgSize, testIm(imgId)); 
%     imgOrg = imnormalize(imgOrg);
    par.imgOrg = imgOrg; par.imgName = imgName;
%     display(['Recover ' par.imgName ' using ' par.recMode]);
    saveFolderText   = ['Result_text' num2str(par.imgSize) '\' ];    if ~exist(saveFolderText, 'dir'); mkdir(saveFolderText);   end;    
    fileNameNave  	 = [saveFolderText imgName '_' par.recMode '_' par.filterMode custorm_note ];
    write_info([fileNameNave  '.txt'], [par.imgName ' size' num2str(par.imgSize)]);
    write_info([fileNameNave  '.txt'], ['               size sigma   iter  sub    time  mse     psnr    ssim   mssim  vsnr   vif    vifp   uqi    ifc    nqm   snr    HFEN    srsim  rfsim  gmsd   fsim   psnrR  psnrG  psnrB psnrAvg psnrSt ssimR  ssimG  ssimB  ssimAv ssim   psnr1 psnr2 psnr3 ssim1 ssim2 ssim3']);

    resutls_all      = cell(1); 
    
    for sub = 1:1:length(par.sparsity)
        subrate = par.sparsity(sub);
        qid_inter = cell(1);
        recImgAll = cell(1);
        
        for trial = 1:1:par.nbrLoop
            display(['Recover ' par.imgName ' using ' par.recMode  ' ' par.filterMode ', subrate' num2str(subrate) ',trial:' num2str(trial) '/' num2str(par.nbrLoop)]);
            
            % -------------- Sensing matrix ------------------------
            [R, G]   = KCS_SensingMtx(par.imgSize, subrate, trial);
            Y        = R*imgOrg*G;
            
            % --------------- Recover the image -----------------------
            tic;           
            opts               = setParams(par.recMode, par);
            [recImg, resInter] = DecWTVNLR(R, G, Y, opts, par.imgOrg, par.recMode);
            recTime(trial)     = toc;
            
            recImgAll{trial} = recImg; 
            
            % save each trial
            qid_inter{trial} = qid_cal(recImg, imgOrg, par.qid_index); 
            write_results([fileNameNave  '.txt'], imgOrg, subrate, 0, qid_inter{trial}, 0, trial, par.nbrLoop);
            
            results.resInter  = resInter;
            results.qid_inter = qid_inter;		
            results.t_org     = recTime;
            results.recImgAll = recImg;
            
            % save the results
            patch3 = ['Results\' par.recMode '_' par.filterMode '_iter' num2str(par.maxIter) '_sub' num2str(subrate*10) ];  if ~exist(patch3, 'dir');  mkdir(patch3); end;
            patch31 = [patch3 '\' par.imgName '_' par.recMode custorm_note '_trial' num2str(trial) '.mat']; save(patch31, 'results', 'par');
        end; % end trial
        % save average
        write_results([fileNameNave  '.txt'], imgOrg, subrate, 0, qid_inter, 0, 0, trial);                
        display(['========== Recovery PSNR:' num2str(qid_inter{trial}.psnr) '============']);
        
		% calculate PSNR of block no boundary				
		results_all{sub} = results;         
        
    end; % end sparsity    
    % save all file 
    save([patch3 '\Full_' par.imgName '_' par.recMode custorm_note '.mat' ], 'results_all', 'par');
end; % end test image
display('END SIMULATION!!!');
