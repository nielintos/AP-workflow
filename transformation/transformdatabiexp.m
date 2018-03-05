function tout=transformdatabiexp(tdata,selch2transf)
tout=tdata;
for ii=1:length(selch2transf)
    dt=tdata{:,selch2transf{ii}};
    [Ti,Wi,Mi,Ai]=estimatebiexponentialparams(dt);
    nn=join(selch2transf{ii});
    disp(['Biexponential params for channels ' nn{1} ':']);
    disp(['  T=' num2str(Ti) '  W=' num2str(Wi) '  M=' num2str(Mi) '  A=' num2str(Ai)]);
    tout{:,selch2transf{ii}}=biexponential(dt,'T',Ti,'W',Wi,'M',Mi,'A',Ai);
end
