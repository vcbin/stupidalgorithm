clc;clear all;
WD='c:\Users\Administrator\Documents\Visual Studio 2010\Projects\PSO\PSO\Input\';
cd(WD);
bst_x_file='bst_x.txt';
before_x_file='ini_pop.txt';
orig_coor_file='original_coor.txt';
node_range=61:300;

bst_x_text=importdata(bst_x_file);
if ( iscell(bst_x_text) )
    bst_x_data=bst_x_text.data;
else
    bst_x_data=bst_x_text;
end
% bst_x_data=load(bst_x_file);
before_x_data=load(before_x_file);
orig_coor_text=importdata(orig_coor_file);
orig_coor_data=orig_coor_text.data;
[m n]=size(orig_coor_data);
item_data=zeros(m,n);
best_before_data=zeros(m,n);
[pop_size dims]=size(before_x_data);

item_node_error=zeros(1,m);
bst_before_node_error=zeros(1,m);
% find best solution within initial population
for i=1:pop_size
    % load one solution
    l=1;
    for j=1:2:dims
        for k=1:n
            item_data(l,k)=before_x_data(i,j+k-1);
        end
        l=l+1;
    end
    % calc avg approx error
    for o=1:m
        item_node_error(o)=norm( item_data(o,:)-orig_coor_data(o,:) );
    end
    avg_err=mean(item_node_error);
    if ( 1==i )
        before_bst_avg_err=avg_err;
        best_before_data=item_data;
        bst_before_node_error=item_node_error;
    else
        if ( avg_err<before_bst_avg_err )
            before_bst_avg_err=avg_err;
            best_before_data=item_data;
            bst_before_node_error=item_node_error;
        end
    end
end

cd('..\Output\');

h=figure;
set(h,'name','Location Comparision','Numbertitle','off');
plot(orig_coor_data(:,1),orig_coor_data(:,2),'go');
hold on;
plot(bst_x_data(:,1),bst_x_data(:,2),'bs');
hold on;
for i=1:m
line([orig_coor_data(i,1) bst_x_data(i,1)],[orig_coor_data(i,2) bst_x_data(i,2)]);
hold on
end
xlabel('X');
ylabel('Y');
legend('original','approximated');
print(gcf,'-dpng','location comparision');

after_node_error=zeros(1,m);
for i=1:m
    after_node_error(i)=norm( bst_x_data(i,:)-orig_coor_data(i,:) );
end
after_avg_err=mean(after_node_error);
% save node error file
fid_before=fopen('node_err_before.txt','w');
for i=1:m
    fprintf(fid_before,'%d\t%f\n',node_range(i),bst_before_node_error(i));
end
fclose(fid_before);
fid_after=fopen('node_err_after.txt','w');
for i=1:m
    fprintf(fid_after,'%d\t%f\n',node_range(i),after_node_error(i));
end
fclose(fid_after);

h=figure;
plot_title='Approximation Error Comparision Per Node';
set(h,'name',plot_title,'Numbertitle','off');
plot(node_range,bst_before_node_error,'g');
% xlabel('#node number');
% ylabel('Approximation Error');
% legend('Before');
% 
% h=figure;
% set(h,'name','Node Approximation Error After','Numbertitle','off');
hold on;
plot(node_range,after_node_error,'b');
xlabel('#node number');
ylabel('Approximation Error');
legend('Before','After');
print(gcf,'-dpng',lower(plot_title));

h=figure;
plot_title='Average Approximation Error Comparision';
set(h,'name',plot_title,'Numbertitle','off');
y=[before_bst_avg_err,after_avg_err];
bar(y);
label_txt={'before','after'};
set(gca,'XTickLabel',label_txt);
for i=1:2
    text(i,y(i),num2str(y(i)),'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center');
end
h = findobj(gca,'Type','patch');
set(h,'EdgeColor','w');
ylabel('Average Approximation Error');
print(gcf,'-dpng',lower(plot_title));