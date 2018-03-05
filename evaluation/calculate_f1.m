
function [t_grp_assign, t_allclustP, t_allclustR, t_allclust] = calculate_f1(t_rs_all,t_clusts_all)
% [t_grp_assign, t_allclustP, t_allclustR, t_allclust] = calculate_f1(t_rs_all,t_clusts_all)
% 
% Compute f1 scores using the Hungarian assignment algorithm. Measurements
% (precision, recall and f1 scores) for each reference standard group 
% (class) and total measurements using all classes.
%
% Outputs:  t_grp_assign (class-cluster assigment together with
%           measurements)
%           t_allclustP (precision for each class and total precision)
%           t_allclustR (recall for each class and total recall)
%           t_allclust (f1 score for each class and total f1 score)
%
% The following parameters need to be specified as arguments:
%
%   Argument     Value
%   ----------   ----------------------------------------------------------
%   t_rs_all     Table with one column representing the reference standard (class#)
%
%   t_clusts_all Table with columns corresponding to different clustering
%                outputs.

t_grp_assign=[];
t_allclustP=[];
t_allclustR=[];
t_allclust=[];

% There should be only one column for RS. Take the first if not
rs=1;
t_rs=t_rs_all(:,rs);
t_grp_truth=grpstats(t_rs,1,'numel');
t_grp_truth(t_grp_truth{:,1}<=0,:)=[];
t_grp_truth.Properties.VariableNames{end}='Truth';

% FOR EACH CLUSTERING SELECTED --> EVALUATE EVERYTHING
t_allclust=t_grp_truth(:,1);
t_allclustP=t_allclust;
t_allclustR=t_allclust;
for cl=1:size(t_clusts_all,2)
    t_cl=t_clusts_all(:,cl);
    t_rscl=[t_rs t_cl];
    % Compute f1 scores
    t_grp=grpstats(t_rscl,1:2);
    t_grp(t_grp{:,1}<=0,:)=[];
    t_grp(t_grp{:,2}<=0,:)=[];
    t_grp.Properties.VariableNames{end}='TP';

    t_grp_detected=grpstats(t_cl,1,'numel');
    t_grp_detected(t_grp_detected{:,1}<=0,:)=[];
    t_grp_detected.Properties.VariableNames{end}='Detected';

    t_grp=join(join(t_grp,t_grp_detected),t_grp_truth);
    t_grp.Precision=t_grp.TP./t_grp.Detected;
    t_grp.Recall=t_grp.TP./t_grp.Truth;
    t_grp.F1=2*(t_grp.Precision.*t_grp.Recall)./(t_grp.Precision+t_grp.Recall);

    % Create matrix to solve the linear assignment problem (minimum cost)
    f1sc=zeros(max(t_grp{:,1}), max(max(t_grp{:,1}),max(t_grp{:,2})));
    for i=1:size(t_grp,1)            
        f1sc(t_grp{i,1},t_grp{i,2})=t_grp.F1(i);
    end
    f1scinv=1-f1sc;

    [rowsol,cost,v,u,costMat] = lapjv(f1scinv);

    idx_assign=ismember(t_grp{:,1:2}, [[1:length(rowsol)]' rowsol(:)],'rows');
    t_grp_assign=t_grp(idx_assign,:);

    % Add classes that were not assigned to a cluster, if any
    t_grp_assign2=outerjoin(t_grp_truth,t_grp_assign,'MergeKeys',1);
    t_grp_assign=t_grp_assign2(:,t_grp_assign.Properties.VariableNames);
    clear t_grp_assign2;
    kk=t_grp_assign{:,3:end};
    kk(isnan(kk))=0;
    t_grp_assign{:,3:end}=kk;

    % Add final row for global numbers (using all classes and all events)
    t_grp_assign{end+1,:}=[Inf,Inf,...
        nansum(t_grp_assign.TP),...
        nansum(t_grp_assign.Detected),...
        nansum(t_grp_assign.Truth),...
        nan,nan,nan];
    t_grp_assign.Precision(end)=t_grp_assign.TP(end)/t_grp_assign.Detected(end);
    t_grp_assign.Recall(end)=t_grp_assign.TP(end)/t_grp_assign.Truth(end);
    t_grp_assign.F1(end)=2*(t_grp_assign.Precision(end)*t_grp_assign.Recall(end))/(t_grp_assign.Precision(end)+t_grp_assign.Recall(end));

    t_allclustP=outerjoin(t_allclust,t_grp_assign(:,[1,end-2]),'MergeKeys',1);
    t_allclustP.Properties.VariableNames{end}=[t_grp_assign.Properties.VariableNames{end-2} '_' t_cl.Properties.VariableNames{1}];

    t_allclustR=outerjoin(t_allclust,t_grp_assign(:,[1,end-1]),'MergeKeys',1);
    t_allclustR.Properties.VariableNames{end}=[t_grp_assign.Properties.VariableNames{end-1} '_' t_cl.Properties.VariableNames{1}];
    
    t_allclust=outerjoin(t_allclust,t_grp_assign(:,[1,end]),'MergeKeys',1);
    t_allclust.Properties.VariableNames{end}=[t_grp_assign.Properties.VariableNames{end} '_' t_cl.Properties.VariableNames{1}];

end

return;          
