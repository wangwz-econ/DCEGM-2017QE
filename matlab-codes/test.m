clearvars
vector = [1;2;3;4;5;5;6;7;7;8;8;9;10;11;12;13;14]
lag = [NaN;vector(1:end-1)];
diff = vector - lag;
if all(diff>=0 | isnan(diff))
  res = vector;
  index = 1:length(vector);
else
  index_start = find(diff<0);
  temp = vector(index_start:end, :);
  res_1 = vector(1:index_start-1);
  res_2 = temp(temp>res_1(end));

  res = [res_1; res_2];
  deleted_index = index_start:1:(index_start + max(find(temp<=res_1(end)) - 1));
  index = setdiff(1:length(vector), deleted_index);
end

[a, b] = fun_refine(vector)
