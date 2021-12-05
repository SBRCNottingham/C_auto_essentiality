function [MCC, TP, TN, FP, FN] = calMCC_mod(trueMat, predictedMat)
% adapted from:
% https://uk.mathworks.com/matlabcentral/fileexchange/47364-true-positives-false-positives-true-negatives-false-negatives-from-2-matrices
% https://itectec.com/matlab/matlab-how-to-define-custom-classification-loss-function/

adder = trueMat + predictedMat
TP = length(find(adder == 2))
TN = length(find(adder == 0))
subtr = trueMat - predictedMat;
FP = length(find(subtr == -1))
FN = length(find(subtr == 1))

if ( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) == 0
    MCC = 0;    % set MCC to zero, if the denominator is zero
else
    MCC = (TP*TN - FP*FN) / ...
        sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) );
end

TP = TP; 
TN = TN; 
FP = FP; 
FN = FN;

end