function Y=biexponential(X, varargin)

    DEFAULT_T=2^(32-1); % T is the formal “top of scale” value. 32 bit storage positive and negative
    DEFAULT_W=0.5;      % W is the number of such decades in the approximately linear region (width parameter). This is a test
    DEFAULT_M=4.5;      % M is the number of decades that the true logarithmic scale approached at the high end of the logicle scale would cover in the plot range. 4.5 is recommended
    DEFAULT_A=0;        % A is the number of additional decades of negative data values to be included. A=0 -->  the most negative value on scale has the same absolute value as the most positive data value in what can be considered the quasilinear region of the scale
    
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
    
    idxnan=find(isnan(X));
    X(idxnan)=0;
    Y=logicleTransform(X,T,W,M,A);
    Y(idxnan)=nan;