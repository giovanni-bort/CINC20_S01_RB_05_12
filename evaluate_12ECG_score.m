% This file contains functions for evaluating algorithms for the 2020 PhysioNet/
% CinC Challenge. You can run it as follows:
%
%   evaluate_12ECG_score(labels, outputs, 'scores.csv')
%
% where 'labels' is a directory containing files with labels, 'outputs' is a
% directory containing files with outputs, and 'scores.csv' (optional) is a
% collection of scores for the outputs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The evaluate_scores function computes a Fbeta measure and a generalizatoin of
% the Jaccard measure but giving missed diagnosis twice as much weight as
% correct diagnoses and false alarms
%
% Inputs:
%   'label_directory' is a directory of comma-delimited text files containing
%   vector of the true labels
%
%   'output_directory' is a directory of comma-delimited text files, where
%   the first row of the file is the predictive labels for each class and
%   the second row of the file is the probability of the class label. 
%   Note that there must be a output for every label.
%
% Outputs:
%
%   'f_beta' is F_beta measure, with beta = 2
%
%   'g_beta' is a generalization of the Jaccard measures but giving missed 
%   diagnoses twice as much weight as correct diagnoses and false alarms, beta = 2
%
%    'accuracy' is accuracy
%
%    'f_measure' is F-measure
%
% Example:
%   Omitted due to length. See the below examples.


function [classes,C_tab,labels, output]= evaluate_12ECG_score(label_directory, output_directory, output_file)

    % Set parameters.
    label_header       = '12ECGLabel';
    output_header  = 'OutputLabel';
    probability_header = 'OutputProbability';

    beta = 2;

    % Find labels and output
    label_files={};

    for f = dir(label_directory)'
        if exist(fullfile(label_directory, f.name), 'file') == 2 && f.name(1) ~= '.' && all(f.name(end - 2 : end) == 'hea')
            label_files{end + 1} = f.name;
        end
    end

    output_files={};
    for f = dir(output_directory)'
        if exist(fullfile(output_directory, f.name), 'file') == 2 && f.name(1) ~= '.' && all(f.name(end - 2 : end) == 'csv')
            output_files{end + 1} = f.name;
        end
    end

    if length(label_files) ~= length(output_files)
        error('Numbers of label and output files must be the same.');
    end


    % Load labels and outputs.
    num_files = length(label_files);
    labels=[];
    output=[];
    output_probabilities=[];

    classes = get_classes(label_directory,label_files);


    for k =1:num_files

	[recording_label,classes_label,single_recording_labels]=get_true_labels([label_directory filesep label_files{k}],classes);

        fid2=fopen(fullfile(output_directory,output_files{k}));
        recording_output = fgetl(fid2);
        classes_output = strsplit(fgetl(fid2),',');
        single_recording_output = cellfun(@str2double,strsplit(fgetl(fid2),','));
     % single_recording_output = [1 1 1 1 1 1 1 1 1]; % TEST CONDITION 
	single_probabilities_output=cellfun(@str2double,strsplit(fgetl(fid2),','));
        fclose(fid2);
if(k<11),fprintf('%4.0f',k,single_recording_output); 
    fprintf('%6.3f',single_probabilities_output);fprintf('\n');end
	% Check labels and outputs for errors.
	if ~(strcmp(classes_label,classes_output))
		error('Numbers of labels and outputs for a file must be the same.');
	end

	if ~(length(single_recording_labels) == length(single_recording_output) || length(single_recording_output) == length(single_probabilities_output))
        	error('Numbers of labels and outputs for a file must be the same.');
	end

	labels = [labels ; single_recording_labels];
	output = [output ; single_recording_output];
	output_probabilities = [output_probabilities ; single_probabilities_output];

    end

    num_classes = length(classes_label);

    K_PROB_MIN=input('if(out=0): min_prob =0.2 , 0.3  2:max_value  3:all_random  4:all_rand_perm_ex[p=0.3] ');
    if(isempty(K_PROB_MIN)),K_PROB_MIN=0.3;end
    K_INI_P=input('K_ini     [1] ');
    if(isempty(K_INI_P)),K_INI_P=1;end
    N_PZ_00=0;
    for IJK=1:num_files
        if(sum(output(IJK,:))==0)
            N_PZ_00=N_PZ_00+1;
            JJJ=find(output_probabilities(IJK,:)>K_PROB_MIN);
           if(K_PROB_MIN==2), [~,JJJ]=max(output_probabilities(IJK,:));end
            output(IJK,JJJ)=1;
        end
        if(K_PROB_MIN==3), 
                JJJ=randperm(9,1);
                output(IJK,:)=0;
                output(IJK,JJJ)=1;
        end
        if(K_PROB_MIN==4), 
                JJJ=randperm_ext(9,1,[12 7 2 9 6 7 18 8 2]);
                output(IJK,:)=0;
                output(IJK,JJJ)=1;
        end
            

       
    end
       fprintf('*** Corrected  %4.0f Output Labels \n',N_PZ_00);
    

    % Compute F_beta measure and the generalization of the Jaccard index
    [accuracy,f_measure,f_beta,g_beta,C_tab] = compute_beta_score(labels(K_INI_P:end,:), output(K_INI_P:end,:), beta, num_classes,classes);


    
    %------------------------------------------------------
    for ijk=1:0 %num_files     % per stampare ECG con nessuna diagnosi
       if(sum(output(ijk,:))==0)
           LAB1=['         '];
           LAB2=['         '];
           ind_LAB1=find(output_probabilities(ijk,:)>0.2);
           ind_LAB2=find(labels(ijk,:)>0);
           LAB1(ind_LAB1)='X';
           LAB2(ind_LAB2)='#';
           fprintf('%4.0f->',ijk);
           for IK=1:9
               fprintf('%6.3f %s%s',output_probabilities(ijk,IK),LAB1(IK),LAB2(IK));end
           fprintf('\n');
       end
    end
    %------------------------------------------------------
    
    % Compute AUC, accuracy, and F-measure.

    [auroc,auprc] = compute_auc(labels, output_probabilities,num_classes);


    % Output results.
    output_string = sprintf('AUROC|AUPRC|Accuracy|F-measure|Fbeta-measure|Gbeta-measure\n%.3f|%.3f|%.3f|%.3f|%.3f|%.3f',...
                             auroc, auprc, accuracy, f_measure, f_beta,g_beta);
                        
    switch nargin
        case 2
            disp(output_string)
        case 3
            fid = fopen(output_file, 'wt');
            fprintf(fid, output_string);
            fclose(fid);
    end


end



% The compute_beta_score function computes the Fbeta-measure giving an specific beta value
% and the G value define at the begining of the file
%
% Inputs:
%   'labels' are the true classes of the recording
%
%   'outputs' are the output classes of your model
%
%   'beta' is the weight
%
% Output:
%   f_beta, Fbeta measure given an specific beta
%
%   g_beta, generalization of the Jaccard measure with a beta weigth

function [accuracy,f_measure,f_beta,g_beta,C_tab] = compute_beta_score(labels, outputs,beta,num_classes,classes)
    % Check inputs for errors.
    % labels : true classifications
    % output:  computed classifications
    % beta:    weight to indices [2]
    % num_classes : number of classes
    % classes: diagnostic classes
    
    fprintf('num_classes=%6.0f\n',num_classes);
    if length(outputs) ~= length(labels)
        error('Numbers of outputs and labels must be the same.');
    end

    [num_recordings,num_classes_from_lab] = size(labels);

        % Check inputs for errors.
    if length(num_classes) ~= length(num_classes_from_lab)
        error('Numbers of classes and labels must be the same.');
    end

    % Populate contingency table.

    fbeta_l = zeros(1,num_classes);
    gbeta_l = zeros(1,num_classes);
    fmeasure_l = zeros(1,num_classes);
    accuracy_l = zeros(1,num_classes);

    f_beta = 0;
    g_beta = 0;
    f_measure = 0;
    accuracy = 0;

    % Weigth function
    C_l = ones(1,num_classes);
TPTN=[];TPTN(num_classes,4)=0;

    for j=1:num_classes
	tp = 0;
	fp = 0;
	fn = 0;
	tn = 0;
	for i = 1 : num_recordings

		num_labels = sum(labels(i,:));

	        if labels(i,j)==1 && outputs(i,j)==1
	            tp = tp + 1/num_labels;
	        elseif labels(i,j)~=1 && outputs(i,j)==1
	            fp = fp + 1/num_labels;
	        elseif labels(i,j)==1 && outputs(i,j)~=1
    		    fn = fn + 1/num_labels;
	        elseif labels(i,j)~=1 && outputs(i,j)~=1
	            tn = tn + 1/num_labels;
	        end
    end

    TPTN(j,1)=tp;
    TPTN(j,2)=fn;
    TPTN(j,3)=fp;
    TPTN(j,4)=tn;
    
	% Summarize contingency table.
        if ((1+beta^2)*tp + (beta*fn) + fp) > 0
	        fbeta_l(j) = ((1+beta^2)*tp) / ((1+beta^2)*tp + (beta^2*fn) + fp);
        else
        	fbeta_l(j) = 1;
        end

	if (tp + (beta*fn) + fp) > 0
	        gbeta_l(j) = tp / (tp + (beta*fn) + fp);
	else
	        gbeta_l(j) = 1;
	end

	if (tp + fp + fn + tn) > 0
	        accuracy_l(j) = (tp+tn) / (tp+fp+fn+tn);
	else
	        accuracy_l(j) = 1;
	end

	if (2*tp + fp + tn) >0
		fmeasure_l(j) = (2*tp)/((2*tp)+fp+fn);
	else
		fmeasure_l(j) = 1;
	end

    end
%     % ----------------------- calcolo e stampa contingency table   ---------------------------
%     C_tab=[];C_tab(1:num_classes,1:num_classes)=0;
%     	for i_p = 1 : num_recordings
%           for j_lab=1:num_classes
%              for j_out=1:num_classes
% 	           if (labels(i_p,j_lab)==1 && outputs(i_p,j_out)==1 ), C_tab(j_lab,j_out)=C_tab(j_lab,j_out)+1;end
%             end
%           end
%         end
%         fprintf('---------------- contingency table old-------------------------------------\n');
%                 fprintf('     -> '); fprintf('%6.0f',1:num_classes);fprintf('\n');
% 
%             for j_lab=1:num_classes,
%                 fprintf('%4.0f -> ',j_lab); fprintf('%6.0f',C_tab(j_lab,:));
%                     fprintf('%6.0f',sum(C_tab(j_lab,:)));   fprintf('\n');
%             end
%              fprintf('     -> ');fprintf('%6.0f',sum(C_tab));fprintf('\n');
%              fprintf('-----------------------------------------------------\n');
%         
%    % --------------------------------------------------------------------------------------------
% ----------------------- calcolo e stampa contingency table   ---------------------------
    C_tab=[];C_tab(1:num_classes,1:num_classes+1)=0;
    	for i_p = 1 : num_recordings
            ZERO=0;
            n_tot_class=sum(labels(i_p,:));
          for j_lab=1:num_classes
	        if (labels(i_p,j_lab)==1 )
                iii=find(outputs(i_p,:)==1);
                if(ismember(j_lab,iii))
                    C_tab(j_lab,j_lab)=C_tab(j_lab,j_lab)+1/n_tot_class;
                else
                    if(numel(iii)>0), C_tab(j_lab,iii)=C_tab(j_lab,iii)+(1/(numel(iii)*n_tot_class));
                    else   C_tab(j_lab,num_classes+1)=C_tab(j_lab,num_classes+1)+1;ZERO=1;
                    end
                end
            end
            
          end
% %         if(ZERO>0),fprintf('zero:%6.0f\n',i_p);end
        end
        fprintf('---------------- contingency table -------------------------------------\n');
                fprintf('      -> '); fprintf('%6.0f',1:num_classes);fprintf('  none   sum   recall precision\n');

            for j_lab=1:num_classes,
                fprintf('%5.0f -> ',j_lab); fprintf('%6.1f',C_tab(j_lab,:));
                    fprintf('%8.0f',sum(C_tab(j_lab,:))); 
                    TP=TPTN(j_lab,1);FN=TPTN(j_lab,2);FP=TPTN(j_lab,3);TN=TPTN(j_lab,4);
                    fprintf('%7.1f%%',100*TP/(TP+FN),100*TP/(TP+FP));
                    
                   fprintf('%8.3f', ((1+beta^2)*TP) / ((1+beta^2)*TP + (beta^2*FN) + FP));
                    fprintf('\n');
            end
             fprintf('  sum -> ');fprintf('%6.0f',sum(C_tab)); fprintf('%8.0f',sum(C_tab(:))); fprintf('\n');
             fprintf('labels-> ');fprintf('%6.0f',sum(labels));fprintf('%14.0f',sum(labels(:)));fprintf('\n');
             fprintf('output-> ');fprintf('%6.0f',sum(outputs));fprintf('%14.0f',sum(outputs(:)));fprintf('\n');
             
             fprintf('-----------------------------------------------------\n');
        
   % --------------------------------------------------------------------------------------------
            
    for i = 1:num_classes
	    f_beta = f_beta + fbeta_l(i)*C_l(i);
            g_beta = g_beta + gbeta_l(i)*C_l(i);
            f_measure = f_measure + fmeasure_l(i)*C_l(i);
            accuracy = accuracy + accuracy_l(i)*C_l(i);
    fprintf('%3.0f %6s ',i,classes{i});          
    fprintf(' F_beta:%8.3f  Acc:%8.3f  G_beta:%8.3f  meas:%8.3f  C:%10.3f \n',fbeta_l(i),accuracy_l(i),gbeta_l(i),fmeasure_l(i),C_l(i));        
    end

    f_beta = f_beta/num_classes;
    g_beta = g_beta/num_classes;
    f_measure = f_measure/num_classes;
    accuracy = accuracy/num_classes;
    
    fprintf('tot summary F_beta:%8.3f  Acc:%8.3f  G_beta:%8.3f  meas:%8.3f \n',f_beta,accuracy,g_beta,f_measure);        
        fprintf('-----------------------------------------------------\n');

end

% The compute_auc function computes AUROC and AUPRC as well as other summary
% statistics (TP, FP, FN, TN, TPR, TNR, PPV, NPV, etc.) that can be exposed
% from this function.
%
% Inputs:
%   'labels' are the true classes of the recording
%
%   'output' are the output classes of your model
%
%   'beta' is the weight
%
%
% Outputs:
%   'auroc' is a scalar that gives the AUROC of the algorithm using its
%   output probabilities, where specificity is interpolated for intermediate
%   sensitivity values.
%
%   'auprc' is a scalar that gives the AUPRC of the algorithm using its
%   output probabilities, where precision is a piecewise constant function of
%   recall.
%

function [auroc, auprc] = compute_auc(labels,probabilities,num_classes)

    % Check inputs for errors.
    if length(probabilities) ~= length(labels)
        error('Numbers of probabilities and labels must be the same.');
    end

    probabilities(isnan(probabilities))=0;


    auroc_l = zeros(1,num_classes);
    auprc_l = zeros(1,num_classes);

    auroc = 0;
    auprc = 0;

    % Weigth function 
    C_l = ones(1,num_classes);

    [num_recordings,num_classes_from_lab] = size(labels);

    for k = 1:num_classes
	    % Find probabilities thresholds.
	    thresholds = flipud(unique(probabilities(:,k)));

	    if thresholds(1) ~= 1
	        thresholds = [1; thresholds];
	    end

	    if thresholds(end) ~= 0
	        thresholds = [thresholds; 0];
	    end

	    m = length(thresholds);

	    % Populate contingency table across probabilities thresholds.
	    tp = zeros(1, m);
	    fp = zeros(1, m);
	    fn = zeros(1, m);
	    tn = zeros(1, m);

	    % Find indices that sort predicted probabilities from largest to smallest.
	    [~, idx] = sort(probabilities(:,k), 'descend');
	

	    i = 1;
	    for j = 1 : m
	        % Initialize contingency table for j-th probabilities threshold.
	        if j == 1
	            tp(j) = 0;
	            fp(j) = 0;
	            fn(j) = sum(labels(:,k));
	            tn(j) = num_recordings - fn(j);
	        else
	            tp(j) = tp(j - 1);
	            fp(j) = fp(j - 1);
	            fn(j) = fn(j - 1);
	            tn(j) = tn(j - 1);
	        end

		% Update contingency table for i-th largest probabilities probability.
	        while i <= num_recordings && probabilities(idx(i),k) >= thresholds(j)
	            if labels(idx(i),k) == 1
	                tp(j) = tp(j) + 1;
	                fn(j) = fn(j) - 1;
	            else
	                fp(j) = fp(j) + 1;
	                tn(j) = tn(j) - 1;
	            end
	            i = i + 1;
	        end
	    end

		

	    % Summarize contingency table.
	    tpr = zeros(1, m);
	    tnr = zeros(1, m);
	    ppv = zeros(1, m);
	    npv = zeros(1, m);

	    for j = 1 : m
	        if tp(j) + fn(j) > 0
	            tpr(j) = tp(j) / (tp(j) + fn(j));
	        else
	            tpr(j) = 1;
	        end

	        if fp(j) + tn(j) > 0
	            tnr(j) = tn(j) / (fp(j) + tn(j));
	        else
	            tnr(j) = 1;
	        end

	        if tp(j) + fp(j) > 0
	            ppv(j) = tp(j) / (tp(j) + fp(j));
	        else
	            ppv(j) = 1;
	        end

	        if fn(j) + tn(j) > 0
	            npv(j) = tn(j) / (fn(j) + tn(j));
	        else
	            npv(j) = 1;
	        end
	    end

	    % Compute AUROC as the area under a piecewise linear function of TPR /
	    % sensitivity (x-axis) and TNR / specificity (y-axis) and AUPRC as the area
	    % under a piecewise constant of TPR / recall (x-axis) and PPV / precision
	    % (y-axis).

	    for j = 1 : m - 1
	        auroc_l(k) = auroc_l(k) + 0.5 * (tpr(j + 1) - tpr(j)) * (tnr(j + 1) + tnr(j));
	        auprc_l(k) = auprc_l(k) + (tpr(j + 1) - tpr(j)) * ppv(j + 1);
	    end
    end


    for i =1:num_classes
	    auroc = auroc + auroc_l(i)*C_l(i);
	    auprc = auprc + auprc_l(i)*C_l(i);
    end

    auroc = auroc/num_classes;
    auprc = auprc/num_classes;
end

% function to obtain the true labels

function [recording_label,classes_label,single_recording_labels]=get_true_labels(input_file,classes)

	classes_label=classes;
	single_recording_labels=zeros(1,length(classes));

	fid=fopen(input_file);
        tline = fgetl(fid);
	tmp_str = strsplit(tline,' ');
	recording_label = tmp_str{1};

        tlines = cell(0,1);
        while ischar(tline)
	        tlines{end+1,1} = tline;
                tline = fgetl(fid);
            if(strncmp(tline,'#Dx',3)>0)  % ***** MODIFIED  VERSION OF startWih *******
  %      	if startsWith(tline,'#Dx')
                        tmp = strsplit(tline,': ');
                        tmp_c = strsplit(tmp{2},',');
                        for j=1:length(tmp_c)
                	        idx2 = find(strcmp(classes,tmp_c{j}));
				single_recording_labels(idx2)=1;
                        end
			break
                end
	end
        fclose(fid);

end


% find unique number of classes
function classes = get_classes(input_directory,files)

        classes={};
        num_files = length(files);
        k=1;
        for i = 1:num_files
                input_file = fullfile(input_directory, files{i});
                fid=fopen(input_file);
                tline = fgetl(fid);
                tlines = cell(0,1);
		while ischar(tline)
                    tlines{end+1,1} = tline;
                    tline = fgetl(fid);
                      if(strncmp(tline,'#Dx',3)>0),  % ***** MODIFIED  VERSION OF startWih *******
                %      if startsWith(tline,'#Dx')
                                tmp = strsplit(tline,': ');
                                tmp_c = strsplit(tmp{2},',');
                                for j=1:length(tmp_c)
                                        idx2 = find(strcmp(classes,tmp_c{j}));
                                        if isempty(idx2)
                                                classes{k}=tmp_c{j};
                                                k=k+1;
                                        end
                                end
                        break
                        end
                end
                fclose(fid);
        end
        classes=sort(classes);
end

