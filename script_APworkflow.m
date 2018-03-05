% Path to R binaries
pthR='C:\Program Files\R\R-3.3.2\bin';

% FCS files to analyze together
% files=[];
files={fullfile(pwd,'Example\SPC_MIX1.fcs'), fullfile(pwd,'Example\SPC_MIX7.fcs')};

% Name of channels to use in the AP workflow  
% chnames=[];
% Example: full-spectral data input mode
chnames={'x488CH1_A','x488CH2_A','x488CH3_A','x488CH4_A','x488CH5_A',...
    'x488CH6_A','x488CH7_A','x488CH8_A','x488CH9_A','x488CH10_A',...
    'x488CH11_A','x488CH12_A','x488CH13_A','x488CH14_A','x488CH15_A',...
    'x488CH16_A','x488CH17_A','x488CH18_A','x488CH19_A','x488CH24_A',...
    'x488CH25_A','x488CH26_A','x488CH27_A','x488CH28_A','x488CH29_A',...
    'x488CH30_A','x488CH31_A','x488CH32_A','x405CH1_A','x405CH2_A',...
    'x640CH1_A','x640CH2_A','x640CH3_A','x640CH4_A','x640CH5_A',...
    'x640CH6_A','x640CH7_A','x640CH8_A','x640CH9_A','x640CH10_A',...
    'x640CH11_A','x640CH12_A','x640CH13_A','x640CH14_A','x640CH15_A',...
    'x640CH16_A','x640CH17_A','x640CH18_A','x640CH19_A','x640CH24_A',...
    'x640CH25_A','x640CH26_A','x640CH27_A','x640CH28_A','x640CH29_A',...
    'x640CH30_A','x640CH31_A','x640CH32_A'};

% Output folder
% savedir=[];
savedir=fullfile(pwd,'Example\X_OUT');

% Type of biexponential transformation: 
%   0 -> no transformation
%   1 -> transform all files together
%   2 -> transform each file individually
% dotransf=[];
dotransf=2;

% Channels to transform with biexponential (logicle)
% namestransfsel2={};
% Example: transformation of all fluorescent channels, 
% each data input mode independently
% 1 -> uncompensated channels
namestransfsel2{1}={'ncV450_A','ncBV711_A','ncPE_A','ncAlexaFluor647_A'};
% 2 -> conventionally compensated channels
namestransfsel2{2}={'ccV450_A','ccBV711_A','ccPE_A','ccAlexaFluor647_A'};
% 3 -> compensated channels by spectral unmixing
namestransfsel2{3}={'scV450_A','scBV711_A','scPE_A','scAlexaFluor647_A'};
% 4 -> raw spectral channels
namestransfsel2{4}={'x488CH1_A','x488CH2_A','x488CH3_A','x488CH4_A','x488CH5_A',...
    'x488CH6_A','x488CH7_A','x488CH8_A','x488CH9_A','x488CH10_A',...
    'x488CH11_A','x488CH12_A','x488CH13_A','x488CH14_A','x488CH15_A',...
    'x488CH16_A','x488CH17_A','x488CH18_A','x488CH19_A','x488CH24_A',...
    'x488CH25_A','x488CH26_A','x488CH27_A','x488CH28_A','x488CH29_A',...
    'x488CH30_A','x488CH31_A','x488CH32_A','x405CH1_A','x405CH2_A',...
    'x640CH1_A','x640CH2_A','x640CH3_A','x640CH4_A','x640CH5_A',...
    'x640CH6_A','x640CH7_A','x640CH8_A','x640CH9_A','x640CH10_A',...
    'x640CH11_A','x640CH12_A','x640CH13_A','x640CH14_A','x640CH15_A',...
    'x640CH16_A','x640CH17_A','x640CH18_A','x640CH19_A','x640CH24_A',...
    'x640CH25_A','x640CH26_A','x640CH27_A','x640CH28_A','x640CH29_A',...
    'x640CH30_A','x640CH31_A','x640CH32_A'};

% Algorithm for subsampling: LocalDensity
% alg=[];
alg='LocalD';

% Number of elements to keep after subsampling (0 if no subsampling)
% n_red=[];
n_red=2500;

% Name of channel with Reference Standard (manual annotations)
% class_str=[];
class_str='CLASS';

% Number of runs
% nTimes=[];
nTimes=2;

path_with_subdirectories = genpath(pwd);
addpath( path_with_subdirectories );

func_APworkflow(files,chnames,savedir,dotransf,namestransfsel2,...
    alg,n_red,class_str,nTimes,pthR);