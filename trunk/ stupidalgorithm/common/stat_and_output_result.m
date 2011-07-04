clc;clear all;
WD='c:\Users\Administrator\Documents\Visual Studio 2010\Projects\PSO\PSO\output\';
%run_stat_file='run_stat_F5_stoptype=stag.txt';
run_stat_file='run_stat.txt';
common_filelist='common_filelist.txt';

para_stat_all_bst_val_filelist='para_stat_all_bst_val_filelist.txt';
para_stat_filelist='para_stat_filelist.txt';

gen_vel_div_filelist='gen_vel_div_filelist.txt';% pso algo specific
mpso_gen_attr_perc_file='mpso_gen_attr_filename.txt'; % mpso algo specific
gen_f_filelist='gen_f_filelist.txt';% de algo specific
gen_pr_filelist='gen_pr_filelist.txt';% de algo specific
gen_eta_filelist='gen_eta_filelist.txt';% ep algo specific

% output file
algo_comp_res='algo_comp_res.txt';% algorithm comparison and ttest result file
para_comp_res='para_comp_res.txt';% parameter combination comparison and ttest result file
seq_algo_desc_res='seq_algo_desc.txt';
de_f_stat_res='de_f_stat.txt';
de_pr_stat_res='de_pr_stat.txt';
ep_eta_stat_res='ep_para_stat.txt';

func_cap={'sphere','griewank','rastrigin','ackey','f5','rosenbrock','step','quartic_with_noise','ws_location','f2','f8','camelback','f12','f13','f3','f4'};
algo_cap={'PSO','mPSO','arPSO','dPSO','*dPSO','PSObc','DE','bbDE','sDE','msDE','spDE','jDE','fEP','ifEP','dgEA'};
line_property=char('b.','go','rx','c+','m*','ys','kd','b^','g<','r>','cp','mh','y-', ...
    'k:','b--','g');
max_color=16;

Min_DE_Algo_Code=7;
Max_DE_Algo_Code=12;
Min_EP_Algo_Code=13;
Max_EP_Algo_Code=14;
DE_Algo_Code=7;
BBDE_Algo_Code=8;
SDE_Algo_Code=9;
MY_SDE_Algo_Code=10;
SPDE_Algo_Code=11;
JDE_Algo_Code=12;
FEP_Algo_Code=13;
IFEP_Algo_Code=14;
Max_Eta_Display_Dim=11;% ep algo specific
ttest_prop=1.0;% t-test threshold proportion

% data column def
func_code_col=1;
algo_code_col=2;
bst_bst_col=3;
wst_bst_col=4;
avg_bst_col=5;
std_col=6;
conv_num_col=7;
conv_perc_col=8;
avg_conv_eval_num_col=9;
avg_time_col=10;
eval_num_col=11;
ob_num_col=12;
avg_stag_gen_col=13;
avg_gen_col=14;

cd (WD);
func_num=0;
func_code_list={};
prev_func_code=-1;
stat=importdata(run_stat_file);
[algo_num data_cols]=size(stat.data);
% find all tested function codes
for i=1:algo_num
    if ( prev_func_code~=stat.data(i,func_code_col) )
        func_code_list{end+1}=stat.data(i,func_code_col);
        prev_func_code=stat.data(i,func_code_col);
        func_num=func_num+1;
    end
end

fid_algo_comp=fopen(algo_comp_res,'w');
fid_para_comp=fopen(para_comp_res,'w');
fid_algo_desc=fopen(seq_algo_desc_res,'w');

dis_item_count=0;% displayed parameter combination counter

fid_de_f_stat=fopen(de_f_stat_res,'w');
fprintf(fid_de_f_stat,'func_code\talgo_code\tf_mean\tf_std\r\n');% print column header
fid_de_pr_stat=fopen(de_pr_stat_res,'w');
fprintf(fid_de_pr_stat,'func_code\talgo_code\tpr_mean\tpr_std\r\n');% print column header
% fprintf(fid_ep_eta_stat,'func_code\talgo_code\teta_mean\teta_std\r\n');% print column header
% [unused func_num]=size(func_code_list);

% excel output vars
row_num=func_num+1;
mean_out=cell(row_num,1);
first_func=true;
writable=true;
prev_algo_num=0;
prev_algo_code_list={};
func_name_list={};
for f_idx=1:func_num % for every tested functon
    cur_func_code=func_code_list{f_idx};
    [func_rows cols]=find( stat.data(:,func_code_col)==cur_func_code );
    cur_data=stat.data(func_rows,:);% get all test data with same test function
    [algo_num data_cols]=size(cur_data);
    % reinitialize stat variables
    one_algo=(algo_num==1);
    has_de=false;
    has_ep=false;
    % generate function name description text
    func_desc=strcat(func_cap{cur_func_code},':');
    func_desc(1)=upper(func_desc(1));
    
    % load common file path list for generational comparison plotting
    i=1;
    cur_idx=1;
    % reinitialize stat variable
    [algo_num unused]=size(cur_data);
    algo_code_list=cell(1,algo_num);
    gen_best_filelist=cell(1,algo_num);
    gen_div_filelist=cell(1,algo_num);
    gen_rad_filelist=cell(1,algo_num);
    all_bst_val_filelist=cell(1,algo_num);
    fidin=fopen(common_filelist,'r');
    while ~feof(fidin) % End Of File?
        tline=fgetl(fidin); % get line from file(no carriage return)
        if ~isempty(tline) % empty line?
            if ( i>1 ) % skip column header
                [cur_func cur_algo_code,cur_all_bst_val_file,cur_gen_best_file,cur_gen_div_file cur_gen_rad_file]= ...
                    strread(tline,'%d\t%d\t%s\t%s\t%s\t%s');
                if ( cur_func==cur_func_code )
                    algo_code_list{cur_idx}=cur_algo_code;
                    all_bst_val_filelist{cur_idx}=cur_all_bst_val_file;
                    gen_best_filelist{cur_idx}=cur_gen_best_file;
                    gen_div_filelist{cur_idx}=cur_gen_div_file;
                    gen_rad_filelist{cur_idx}=cur_gen_rad_file;
                    cur_idx=cur_idx+1;
                end
            end
            i=i+1;
        end
    end
    fclose(fidin);
    
    if ( first_func )
        prev_algo_num=algo_num;
        prev_algo_code_list=algo_code_list;
    end
    
    % generate X axis labels
    x_algo_cap_orig=cell(1,algo_num);
    x_algo_cap=x_algo_cap_orig;
    de_algo_num=0;
    ep_algo_num=0;
    same_algo_idx=1;
    first_same=true;
    for i=1:algo_num
        x_algo_cap_orig(i)=algo_cap(cur_data(i,algo_code_col));
        x_algo_cap(i)=x_algo_cap_orig(i);
        % label different parameter combination of SAME algorithm with
        % format: algo_name(1),algo_name(2),
        if ( i>1 )
            if ( strcmp(x_algo_cap_orig(i),x_algo_cap_orig(i-1)) )
                if ( first_same )
                    x_algo_cap{i-1}=strcat(x_algo_cap_orig{i-1},'(',num2str(same_algo_idx),')');
                    x_algo_cap{i}=strcat(x_algo_cap_orig{i},'(',num2str(same_algo_idx+1),')');
                    same_algo_idx=same_algo_idx+2;
                    first_same=false;
                else
                    x_algo_cap{i}=strcat(x_algo_cap_orig{i},'(',num2str(same_algo_idx),')');
                    same_algo_idx=same_algo_idx+1;
                end
            else
                same_algo_idx=1;
                first_same=true;
            end
        end
        if ( cur_data(i,algo_code_col)>=Min_DE_Algo_Code && cur_data(i,algo_code_col)<=Max_DE_Algo_Code )
            has_de=true;
            de_algo_num=de_algo_num+1;
        end
        if ( cur_data(i,algo_code_col)>=Min_EP_Algo_Code && cur_data(i,algo_code_col)<=Max_EP_Algo_Code )
            has_ep=true;
            if ( cur_data(i,algo_code_col)==FEP_Algo_Code )
                has_fep=true;
            end
            if ( cur_data(i,algo_code_col)==IFEP_Algo_Code )
                has_ifep=true;
            end
            ep_algo_num=ep_algo_num+1;
        end
    end % for i=1:algo_num
    
    % record best algo number for performance comparison and analysis
    avg_bst=cur_data(:,avg_bst_col);
    [best_val best_idx]=min(avg_bst);
    best_algo_code=cur_data(best_idx,algo_code_col);
    equal_count=0;% ttest counter
    equ_less_desc={};
    equ_more_desc={};
    equ_less_num_list={};
    equ_more_num_list={};
    equ_P_list={};
    best_equ_count=0;
    best_count=1;
    best_less_desc={};
    best_more_desc={};
    best_more_avg={};
    best_P_list={};
    max_val=0;
    min_val=0;
    has_no_optimum=false;% indicator of optimum algorithm within algorithm portfolio's existence
    has_same_algo=false;
    bst_desc_txt='';
    
    % set has_sequential_value flag for every algorithm instance
    has_seq_val=zeros(1,algo_num);
    seq_desc_idx=cell(1,algo_num);% mapping of algorithm code to description text index
    desc_count=1;
    first_same=true;
    [item_count unused]=size(cur_data);
    for i=1:item_count-1
        if ( cur_data(i,algo_code_col)==cur_data(i+1,algo_code_col) )
            has_seq_val(i)=true;
            has_seq_val(i+1)=true;
            seq_desc_idx{i}=desc_count;
            seq_desc_idx{i+1}=desc_count+1;
            desc_count=desc_count+2;
        else
            if ( i==1 )
                has_seq_val(1)=false;
                seq_desc_idx{1}=-1;
            end
            has_seq_val(i+1)=false;
            seq_desc_idx{i+1}=-1;
        end
    end
    
    % load parameter description text
    line_count=1;
    para_desc_text=cell(1,algo_num);
    fid_para=fopen(para_stat_all_bst_val_filelist);
    while ~feof(fid_para) % End Of File?
        tline=fgetl(fid_para); % get line from file(no carriage return)
        if ~isempty(tline) && ischar(tline) % empty line?
            [cur_func cur_para_desc unused]=strread(tline,'%d\t%s\t%s');
            if ( cur_func==cur_func_code )
                para_desc_text{line_count}=cur_para_desc;
                line_count=line_count+1;
            end
        end
    end
    fclose(fid_para);
    
    same_algo_idx=1;
    first_same=true;
    for i=1:algo_num
        x_algo_cap_orig(i)=algo_cap(cur_data(i,algo_code_col));
        x_algo_cap(i)=x_algo_cap_orig(i);
        % label different parameter combination of SAME algorithm with
        % format: algo_name(1),algo_name(2),
        if ( i>1 )
            if ( strcmp(x_algo_cap_orig(i),x_algo_cap_orig(i-1)) )
                if ( first_same )
                    x_algo_cap{i-1}=strcat(x_algo_cap_orig{i-1},'(',num2str(same_algo_idx),')');
                    x_algo_cap{i}=strcat(x_algo_cap_orig{i},'(',num2str(same_algo_idx+1),')');
                    same_algo_idx=same_algo_idx+2;
                    first_same=false;
                else
                    x_algo_cap{i}=strcat(x_algo_cap_orig{i},'(',num2str(same_algo_idx),')');
                    same_algo_idx=same_algo_idx+1;
                end
            else
                same_algo_idx=1;
                first_same=true;
            end
        end
    end % for i=1:algo_num
    
    % output algo_name(1,...,n) -> parameter value description for
    % plot reading convenience
    [unused cols]=find( has_seq_val(:)==true );
    if ( false==isempty(cols) )
        func_name=func_cap{cur_func_code};
        fprintf(fid_algo_desc,'Function %s\n',func_name);
        same_algo_idx=1;
        first_same=true;
        for i=1:algo_num
            if ( true==has_seq_val(i) )
                cur_algo_code=cur_data(i,algo_code_col);
                % find corresponding entry index in para desc file
                [rows unused]=find ( cur_data(1:i,algo_code_col)==cur_algo_code );
                tmp_data=cur_data(rows,:);
                [dst_idx unused]=size(tmp_data);
                fprintf(fid_algo_desc,'%s:%s\r\n',x_algo_cap{i},char(para_desc_text{dst_idx}));
            end
        end % for i=1:algo_num
        fprintf(fid_algo_desc,'\r\n');
    end
    
    % record best choice parameter description
    if ( true==has_seq_val(best_idx) )
        % find corresponding entry index in para desc file
        [rows unused]=find ( cur_data(1:best_idx,algo_code_col)==best_algo_code );
        tmp_data=cur_data(rows,:);
        [dst_idx unused]=size(tmp_data);
        best_para_desc_text=para_desc_text{dst_idx};
    end
    
    % start mean gbest comparison
    for i=1:algo_num-1
        for j=i+1:algo_num
            min_val=avg_bst(i);
            max_val=avg_bst(j);
            less_idx=i;
            more_idx=j;
            if ( avg_bst(j)<avg_bst(i) )
                min_val=avg_bst(j);
                max_val=avg_bst(i);
                less_idx=j;
                more_idx=i;
            end
            
            if ( (max_val-min_val) <= ttest_prop*abs(min_val) ) % perform t-test if mean gbest value are close(<20%)
                for k=1:algo_num
                    % find corresponding data filename
                    if ( k~=less_idx && k~=more_idx )
                        continue;
                    end
                    file_path=char(all_bst_val_filelist{k});
                    all_bst_val_data=load(file_path);
                    if ( k==less_idx )
                        less_all_bst_vals=all_bst_val_data;
                    else
                        more_all_bst_vals=all_bst_val_data;
                    end
                end % for k=1:algo_num
                [H P CI]=ttest2(less_all_bst_vals,more_all_bst_vals,0.05,'both','unequal');
                identical=all(less_all_bst_vals==more_all_bst_vals);
                if ( 0==H || identical )
                    equal_count=equal_count+1;
                    % find corresponding entry index in para desc file for
                    % less indx
                    [rows unused]=find ( cur_data(1:less_idx,algo_code_col)==algo_code_list{less_idx} );
                    tmp_data=cur_data(rows,:);
                    [dst_idx unused]=size(tmp_data);
                    if ( has_seq_val(less_idx) )
                        less_desc_text=strcat(algo_cap(algo_code_list{less_idx}),'(',para_desc_text{dst_idx},')');
                    else
                        less_desc_text=algo_cap(algo_code_list{less_idx});
                    end
                    % find corresponding entry index in para desc file for
                    % more_idx
                    [rows unused]=find ( cur_data(1:more_idx,algo_code_col)==algo_code_list{more_idx} );
                    tmp_data=cur_data(rows,:);
                    [dst_idx unused]=size(tmp_data);
                    if ( has_seq_val(more_idx) )
                        more_desc_text=strcat(algo_cap(algo_code_list{more_idx}),'(',para_desc_text{dst_idx},')');
                    else
                        more_desc_text=algo_cap(algo_code_list{more_idx});
                    end
                    equ_less_desc{end+1}=less_desc_text;
                    equ_more_desc{end+1}=more_desc_text;
                    if ( 0==H )
                        t_P=P;
                    else % two group data are identical
                        t_P=1.0;
                    end
                    equ_P_list{end+1}=t_P;
                    if ( less_idx==best_idx )
                        best_equ_count=best_equ_count+1;
                        best_count=best_count+1;
                        has_no_optimum=true;
                        best_less_desc{end+1}=less_desc_text;
                        best_more_desc{end+1}=more_desc_text;
                        best_P_list{end+1}=t_P;
                        best_more_avg{end+1}=cur_data(more_idx,avg_bst_col);
                    end
                end % if ( 0==H )
            end %if ( (max-min)>ttest_prop*min )
        end% for j=i+1:algo_num
    end% for i=1:algo_num-1
    
    % output algorithm choosing result
    algo_name=char(algo_cap(best_algo_code));
    if ( has_seq_val(best_idx)==true )
        best_desc_text=sprintf('%s of %s',char(best_para_desc_text),algo_name);
    else
        best_desc_text=sprintf('%s',algo_name);
    end
    % output single best algorithm choice
    func_name=func_cap{cur_func_code};
    func_name(1)=upper(func_name(1));
    if ( best_count==1 )
        fprintf(fid_algo_comp,'%d. %s is the optimum algorithm among all tested algorithms for function %s,mean of gbest=%f.\r\n',...
            f_idx,best_desc_text, ...
            func_name,best_val);
    else
        fprintf(fid_algo_comp,'%d. NO single optimum algorithm among all tested algorithms for function %s\r\n',...
            f_idx,func_name);
    end
    % output multiple best algorithm choices if exist
    if ( best_equ_count>0 )
        fprintf(fid_algo_comp,'\tbest algorithm:\r\n');
        best_info=struct('best_less_desc',best_less_desc,'best_more_desc',best_more_desc,'P',best_P_list,'best_more_avg',best_more_avg);
        diff_text='two mean''s difference are statistical INsignificant';
        [unused order]=sort([best_info.best_more_avg],'ascend');
        best_info=best_info(order);
        for i=1:best_equ_count
            fprintf(fid_algo_comp,'\t\t(%d).\t%s(mean=%f) <-> %s(mean=%f): %s,P=%.1f%%\r\n',i,char(best_info(i).best_less_desc),best_val,...
                char(best_info(i).best_more_desc),best_info(i).best_more_avg, ...
                diff_text, ...
                best_info(i).P*100);
        end
        best_ratio=((best_equ_count+1)/algo_num)*100;
        fprintf(fid_algo_comp,'\ttotal algorithm number:%d,best_ratio=%.2f%%\r\n',algo_num,best_ratio);
        fprintf(fid_algo_comp,'\r\n');
    end % if ( best_count>0 )
    % output ttest comparison result
    if ( equal_count>0 )
        fprintf(fid_algo_comp,'\tequal algorithm:\r\n');
        equal_info=struct('equ_less_desc',equ_less_desc,'equ_more_desc',equ_more_desc,'P',equ_P_list);
        diff_text='two mean''s difference are statistical INsignificant';
        [unused order]=sort([equal_info.P],'descend');
        equal_info=equal_info(order);
        for i=1:equal_count
            fprintf(fid_algo_comp,'\t\t(%d).\t%s <-> %s: %s,P=%.1f%%\r\n',i,char(equal_info(i).equ_less_desc),...
                char(equal_info(i).equ_more_desc), ...
                diff_text, ...
                equal_info(i).P*100);
        end
        rel_num=(algo_num*(algo_num-1))/2;
        if ( rel_num > 0 )
            equal_ratio=(equal_count/rel_num)*100;
        else
            equal_ratio=0.0;
        end
        fprintf(fid_algo_comp,'\ttotal relationship number:%d,equal_ratio=%.2f%%\r\n\r\n',rel_num,equal_ratio);
    end % if ( equal_count>0 )
    
    % start sequential parameter value combination analysis
    dis_item_count=0;
    fid_para=fopen(para_stat_filelist,'r');
    while ~feof(fid_para) % End Of File?
        tline=fgetl(fid_para); % get line from file(no carriage return)
        if ( ~isempty(tline) && ischar(tline) ) % empty line?
            equal_count=0;% reinitialize ttest equal counter
            has_no_optimum=false;
            equ_less_num_list={};
            equ_more_num_list={};
            equ_less_desc={};
            equ_more_desc={};
            equ_P_list={};
            best_equ_count=0;
            best_count=1;
            best_less_num_list={};
            best_more_num_list={};
            best_less_desc={};
            best_more_desc={};
            best_more_avg={};
            best_P_list={};
            
            best_desc_txt='';
            dis_item_count=dis_item_count+1;
            [cur_func algo_type,para_stat_file_path]=strread(tline,'%d\t%d\t%s');
            if ( cur_func==cur_func_code )
                para_stat_file_path=char(para_stat_file_path);
                para_stat=importdata(para_stat_file_path);
                % parsing file to get number of sequential variables
                seq_var_num_txt=char(para_stat.textdata(2,:));
                pos=regexp(seq_var_num_txt,'[0-9]+');
                seq_var_num_txt=seq_var_num_txt(pos:end);
                seq_var_num=str2num(seq_var_num_txt);
                stat_data=para_stat.data;
                % get mean gbest column data for performance comparison
                avg_bst_offset=3;
                avg_bst_col=seq_var_num+avg_bst_offset;
                avg_bst=para_stat.data(:,avg_bst_col);
                [best_val best_idx]=min(avg_bst);
                % split sequential variable name
                seq_txt=char(para_stat.textdata(1,:));
                var_name=regexp(seq_txt,'\t','split');
                var_name=var_name(1:seq_var_num);
                [comb_num n]=size(para_stat.data);
                seq_var_val=zeros(1,seq_var_num);
                % start mean gbest comparison of different parameter combination
                for i=1:comb_num-1
                    % assemble 'less' combination output text
                    cur_seq_out='';
                    for j=1:seq_var_num
                        cur_val=para_stat.data(i,j);
                        cur_seq_out=strcat(cur_seq_out,var_name(j),'=',num2str(cur_val));
                        if ( j~=seq_var_num )
                            cur_seq_out=strcat(cur_seq_out,',');
                        end
                    end
                    
                    for j=i+1:comb_num
                        min_val=avg_bst(i);
                        max_val=avg_bst(j);
                        less_idx=i;
                        more_idx=j;
                        less_seq_out=cur_seq_out;
                        more_seq_out='';
                        % assemble 'more' combination output text
                        for k=1:seq_var_num
                            cur_val=para_stat.data(j,k);
                            more_seq_out=strcat(more_seq_out,var_name(k),'=',num2str(cur_val));
                            if ( k~=seq_var_num )
                                more_seq_out=strcat(more_seq_out,',');
                            end
                        end
                        
                        if ( max_val<min_val )
                            min_val=avg_bst(j);
                            max_val=avg_bst(i);
                            less_idx=j;
                            more_idx=i;
                            % swap output text for outputing order
                            tmp=less_seq_out;
                            less_seq_out=more_seq_out;
                            more_seq_out=tmp;
                        end
                        % record best value combination description
                        if ( less_idx==best_idx )
                            best_desc_txt=less_seq_out;
                        end
                        
                        % load ttest sample data accordingly
                        fidin=fopen(para_stat_all_bst_val_filelist,'r');
                        p=1;% line counter
                        while ~feof(fidin) % End Of File?
                            tline=fgetl(fidin); % get line from file(no carriage return)
                            if ~isempty(tline) % empty line?
                                [func_code para_text all_bst_val_file_path]=strread(tline,'%d\t%s\t%s');
                                if ( func_code==cur_func_code )
                                    if ( strcmp(para_text,less_seq_out) )
                                        file_path=char(all_bst_val_file_path);
                                        all_bst_val_data=load(file_path);
                                        less_all_bst_vals=all_bst_val_data;
                                    end
                                    if ( strcmp(para_text,more_seq_out) )
                                        file_path=char(all_bst_val_file_path);
                                        all_bst_val_data=load(file_path);
                                        more_all_bst_vals=all_bst_val_data;
                                    end
                                    p=p+1;
                                end % if ( func_code==cur_func_code )
                            end % if ~isempty(tline)
                        end % while ~feof(fidin)
                        fclose(fidin);
                        
                        if ( (max_val-min_val) <= ttest_prop*abs(min_val) )
                            [H P CI]=ttest2(less_all_bst_vals,more_all_bst_vals,0.05,'both','equal');
                            identical=all(less_all_bst_vals==more_all_bst_vals);
                            if ( 0==H || identical )
                                % record equal case information
                                equal_count=equal_count+1;
                                equ_less_num_list{end+1}=less_idx;
                                equ_more_num_list{end+1}=more_idx;
                                equ_less_desc{end+1}=char(less_seq_out);
                                equ_more_desc{end+1}=char(more_seq_out);
                                if ( 0==H )
                                    t_P=P;
                                else
                                    t_P=1.0;
                                end
                                equ_P_list{end+1}=t_P;
                                if ( less_idx==best_idx )
                                    best_equ_count=best_equ_count+1;
                                    best_count=best_count+1;
                                    has_no_optimum=true;
                                    best_less_num_list{end+1}=less_idx;
                                    best_more_num_list{end+1}=more_idx;
                                    best_less_desc{end+1}=char(less_seq_out);
                                    best_more_desc{end+1}=char(more_seq_out);
                                    best_more_avg{end+1}=avg_bst(more_idx);
                                    best_P_list{end+1}=t_P;
                                end
                            end % if ( 0==H )
                        end %if ( (max-min)>0.1*min )
                    end % for j=i+1:comb_num
                end % for i=1:comb_num-1
                
                % output single optimum parameters combination
                algo_name=char(algo_cap(algo_type));
                func_name=func_cap{cur_func_code};
                func_name(1)=upper(func_name(1));
                if ( best_count==1 )
                    fprintf(fid_para_comp,'%d. %s is the optimum combination for %s algorithm for function %s,mean of gbest=%f.\r\n',...
                        dis_item_count,char(best_desc_txt),algo_name,func_name,best_val);
                else
                    fprintf(fid_para_comp,'%d. NO single optimum combination for %s algorithm for function %s\r\n',dis_item_count, ...
                        algo_name,func_name);
                end
                % output multiple best combination set if it exist
                if ( best_equ_count>0 )
                    fprintf(fid_para_comp,'\tbest parameter combination:\r\n');
                    best_info=struct('best_less_desc',best_less_desc,'best_more_desc',best_more_desc,'best_less_num_list',best_less_num_list, ...
                        'best_more_num_list',best_more_num_list,'P',best_P_list,'best_more_avg',best_more_avg);
                    diff_text='two mean''s difference are statistical INsignificant';
                    [unused order]=sort([best_info.best_more_avg],'ascend');
                    best_info=best_info(order);
                    for i=1:best_equ_count
                        fprintf(fid_para_comp,'\t\t(%d).\t%s(mean=%f) <-> %s(mean=%f): %s,P=%.1f%%\r\n',i,best_info(i).best_less_desc,best_val, ...
                            best_info(i).best_more_desc,best_info(i).best_more_avg, ...
                            diff_text, ...
                            best_info(i).P*100);
                    end
                    best_ratio=((best_equ_count+1)/comb_num)*100;
                    fprintf(fid_para_comp,'\ttotal parameter combination number:%d,best_ratio=%.2f%%\r\n',comb_num,best_ratio);
                    fprintf(fid_para_comp,'\r\n');
                end % if ( best_count>0 )
                % output ttest comparison result
                if ( equal_count>0 )
                    fprintf(fid_para_comp,'\tequal parameter combination:\r\n');
                    equal_info=struct('equ_less_desc',equ_less_desc,'equ_more_desc',equ_more_desc,'equ_less_num_list',equ_less_num_list, ...
                        'equ_more_num_list',equ_more_num_list,'P',equ_P_list);
                    diff_text='two mean''s difference are statistical INsignificant';
                    [unused order]=sort([equal_info.P],'descend');
                    equal_info=equal_info(order);
                    for i=1:equal_count
                        fprintf(fid_para_comp,'\t\t(%d).\t%s <-> %s: %s,P=%.1f%%\r\n',i,equal_info(i).equ_less_desc,equal_info(i).equ_more_desc, ...
                            diff_text, ...
                            equal_info(i).P*100);
                    end
                    rel_num=(comb_num*(comb_num-1))/2;
                    if ( rel_num > 0 )
                        equal_ratio=(equal_count/rel_num)*100;
                    else
                        equal_ratio=0.0;
                    end
                    fprintf(fid_para_comp,'\ttotal relationship number:%d,equal_ratio=%.2f%%\r\n\r\n',rel_num,equal_ratio);
                end % if ( equal_count>0 )
            end % if ( cur_func==cur_func_code )
        end  % if ~isempty(tline) % empty line?
    end % while ~feof(fid_para)
    fclose(fid_para);
    
    % output DE specific parameters stat value to files
    if ( has_de )
        fidin=fopen(gen_f_filelist,'r');
        while ~feof(fidin) % End Of File?
            tline=fgetl(fidin); % get line from file(no carriage return)
            if ~isempty(tline) % empty line?
                [cur_func algo_type,gen_f_file_path]=strread(tline,'%d\t%d\t%s');
                if ( cur_func==cur_func_code  )
                    gen_f_file_path=char(gen_f_file_path);
                    gen_f_data=load(gen_f_file_path);
                    fprintf(fid_de_f_stat,'%d\t%d\t%f\t%f\r\n',cur_func_code,algo_type,mean(gen_f_data(:,2)),std(gen_f_data(:,2)));% output f value stat to file
                end %if ( cur_func==cur_func_code  )
            end% if ~isempty(tline)
        end
        fclose(fidin);
        
        fidin=fopen(gen_pr_filelist,'r');
        while ~feof(fidin) % End Of File?
            tline=fgetl(fidin); % get line from file(no carriage return)
            if ~isempty(tline) % empty line?
                [cur_func algo_type,gen_pr_file_path]=strread(tline,'%d\t%d\t%s');
                if ( cur_func==cur_func_code )
                    gen_pr_file_path=char(gen_pr_file_path);
                    gen_pr_data=load(gen_pr_file_path);
                    fprintf(fid_de_pr_stat,'%d\t%d\t%f\t%f\r\n',cur_func_code,algo_type,mean(gen_pr_data(:,2)),std(gen_pr_data(:,2)));% output pr value stat to file
                end %if ( cur_func==cur_func_code  )
            end% if ~isempty(tline)
        end
        fclose(fidin);
    end % if ( has_de )
    
    % start excel output writing
    cur_row=f_idx+1;
    if ( first_func )
        % print column header for excel output
        col_num=algo_num*2+1;
        mean_out{1}=cell(1,col_num);
        mean_out{1,1}='';
        for i=1:algo_num
            mean_out{1,i*2}=[x_algo_cap{i},' ','mean'];
            mean_out{1,i*2+1}=[x_algo_cap{i},' ','std'];
        end
        conv_out=mean_out;
        first_func=false;
    end
    % check colomn number equality
    if ( algo_num~=prev_algo_num )
        writable=false;
    else
        for i=1:algo_num
            if ( algo_code_list{i}~=prev_algo_code_list{i} )
                writable=false;
                break;
            end
        end
    end
    if ( writable )
        % collect algorithm performance indices for excel output
        [ unused func_name_len]=size(func_desc);
        func_name=func_desc(1:func_name_len-1);% delete ':'
        mean_out{cur_row}=cell(1,col_num);
        mean_out{cur_row,1}=func_name;
        conv_out{cur_row}=mean_out{cur_row};
        func_name_list{end+1}=func_name;
        for i=1:algo_num
            tmp_data=num2str(cur_data(i,avg_bst_col));
            mean_out{cur_row,i*2}=tmp_data;
            mean_out{cur_row,i*2+1}=num2str(cur_data(i,std_col));
            tmp_data=strcat(num2str(cur_data(i,conv_num_col)),'(',num2str(cur_data(i,conv_perc_col)*100),'%)');
            conv_out{cur_row,i+1}=tmp_data;
        end
    end
    prev_algo_code_list=algo_code_list;
    prev_algo_num=algo_num;
end %for f_idx=1:func_num
fclose(fid_algo_desc);
fclose(fid_para_comp);
fclose(fid_algo_comp);
fclose(fid_de_f_stat);
fclose(fid_de_pr_stat);
% de_f_stat,de_pr_stat files done.

if ( writable )
    xlswrite('Results.xls',mean_out,'mean&std');
    xlswrite('Results.xls',conv_out,'convergence');
    if ( has_de )
        f_out=cell(row_num,col_num);
        f_out(1,:)=mean_out(1,:);
        f_out(:,1)=mean_out(:,1);
        pr_out=f_out;
        f_mean_col=3;
        f_std_col=4;
        pr_mean_col=3;
        pr_std_col=4;
        f_res=importdata(de_f_stat_res);
        pr_res=importdata(de_pr_stat_res);
        % another loop for f,pr excel output
        for i=1:func_num % for every tested functon
            cur_func_code=func_code_list{i};
            [func_rows cols]=find( f_res.data(:,func_code_col)==cur_func_code );
            cur_f_data=f_res.data(func_rows,:);% get all test data with same test function
            [func_rows cols]=find( pr_res.data(:,func_code_col)==cur_func_code );
            cur_pr_data=pr_res.data(func_rows,:);% get all test data with same test function
            for j=1:algo_num
                f_out{i+1,j*2}=num2str(cur_f_data(j,f_mean_col));
                f_out{i+1,j*2+1}=num2str(cur_f_data(j,f_std_col));
                pr_out{i+1,j*2}=num2str(cur_pr_data(j,pr_mean_col));
                pr_out{i+1,j*2+1}=num2str(cur_pr_data(j,pr_std_col));
            end
        end
        xlswrite('Results.xls',f_out,'f value');
        xlswrite('Results.xls',pr_out,'pr value');
    end
end % if ( writable )
