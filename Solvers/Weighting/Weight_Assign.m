function [Wx, Wy, rec_filtered] = Weight_Assign(weight_mode_id, rec, rec_dx, rec_dy, tau1, sigma)
% Function to calculate the weights for total variaiotni
%   [Wx, Wy] = Weight_Assign(weight_mode, rec_spatiall, rec_dx, rec_dy, tau1, tau2)
%   - Input:
%       + weight_mode: NO, ENOW, HENOW, EDGE-CS, CANDES
%       + rec_dx     : input gradient image dx
%       + rec_dy     : input gradient image dy
%       + tau1       : threshold 1; which is used in ENOW, HENOW
%       + tau2       : threshold 2; which is used in EdgeCS
weight_mode     = {'NO', 'ENOW', 'HENOW', 'edgeCS', 'Candes' };
[dx, dy] = gradCal3(rec, 3);
switch weight_mode{weight_mode_id}
    case 'NO'
        Wx      = ones(size(rec_dx));
        Wy      = Wx; 
        
end;
rec_filtered  = rec;