% ESTIMATION EXTRACTED FROM CYTOFKIT IN R
% https://github.com/Bioconductor-mirror/cytofkit/blob/master/R/cytof_preProcess.R
function [T,W,M,A]=estimatebiexponentialparams_v3(X, varargin)   
    X=X(:);
    DEFAULT_T=quantile(X(:),1-1e-4);%nanmax(X(:)); % T is the formal “top of scale” value. 32 bit storage positive and negative
    %DEFAULT_W=0.5;      % W is the number of such decades in the approximately linear region (width parameter). This is a test
    DEFAULT_M=4.5;      % M is the number of decades that the true logarithmic scale approached at the high end of the logicle scale would cover in the plot range. 4.5 is recommended
    DEFAULT_A=0;        % A is the number of additional decades of negative data values to be included. A=0 -->  the most negative value on scale has the same absolute value as the most positive data value in what can be considered the quasilinear region of the scale
    idx_neg=find(X<0);
    q1=quantile(X(idx_neg), 0.25);
    q3=quantile(X(idx_neg), 0.75);
    iqr=q3-q1;
    idx_neg2=idx_neg(X(idx_neg)>=(q1-1.5*iqr));
    if isempty(idx_neg2)
        DEFAULT_W=0;
    else
        r=quantile(X(idx_neg2),0.05)+eps;
        if (10^DEFAULT_M*abs(r)) <= DEFAULT_T
            DEFAULT_W=0;
        else
            DEFAULT_W=(DEFAULT_M-log10(DEFAULT_T/abs(r)))/2;
            if DEFAULT_W>2
                DEFAULT_W=2;
            end
        end
    end
    
    p = inputParser;
    addRequired(p,'X',@isnumeric);
    addParamValue(p,'T',DEFAULT_T,@isnumeric);
    addParamValue(p,'W',DEFAULT_W,@isnumeric);
    addParamValue(p,'M',DEFAULT_M,@isnumeric);
    addParamValue(p,'A',DEFAULT_A,@isnumeric);

    parse(p,X, varargin{:});

    W = p.Results.W;
    M = p.Results.M;
    T = p.Results.T;
    A = p.Results.A;