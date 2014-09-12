function evaluate_results()

results = load('temp_res');

classes = (unique(results(:,1)));

sum_acc = 0;
for i = classes'
	cur_class_acc = mean(results(find(results(:,1)==i),2) == i)
	sum_acc = sum_acc + cur_class_acc;
end

uar = sum_acc/length(classes)
