function two_sample_t_test(data_1,data_2)

% Unpaired t-test
[h, p, ci, stats] = ttest2(data_1, data_2);

% Bar graph
figure;
means = [mean(data_1), mean(data_2)];
errors = [std(data_1)/sqrt(length(data_1)), std(data_2)/sqrt(length(data_2))];
bar(means);
hold on;
errorbar(means, errors, '.');
hold off;

% Display p-value
p_value_text = sprintf('p = %.3f', p);
text(1.5, max(means) + max(errors), p_value_text, 'HorizontalAlignment', 'center');

% Labels and Title
xlabel('Group');
ylabel('Mean Value');
title('Unpaired t-test and Bar Graph with P-value');

end