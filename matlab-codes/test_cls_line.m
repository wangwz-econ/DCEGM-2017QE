%%
clear all

x1 = [1, 2, 3, 4, 5, 6, 7];
y1 = x1 - 2;

x2 = [2;3;4;5;6;10];
y2 = x2 * 2;

x3 = linspace(0,1,100);
y3 = log(x3);

x4 = linspace(1, 2, 1000)';
y4 = exp(x4);

obj1 = cls_line(x1, y1);
obj2 = cls_line(x2, y2);
obj3 = cls_line(x3, y3);
obj4 = cls_line(x4, y4);
obj_array1 = [obj1, obj2];
% obj_array1.cls_len

obj_array2 = [obj1, obj2; obj3, obj4];
% obj_array2.cls_len

obj_array3 = [obj1; obj2; obj3];
% obj_array3.cls_len

obj_array4 = repmat([obj1, obj2], 1, 2);
size(obj_array4)


union(obj_array1.x)
union(obj_array2.x)
union(obj_array3.x)
obj_array4.x
%%
a = cell(2,1)
a{1} = [1:4]
a{2} = [2:5]
cell2mat(a)

%%

ojb1.cls_grow(obj3, true)
cls_grow(obj1, obj3, true)

obj_array_grow1 = [obj1, obj2];
obj_array_grow2 = [obj3, obj4];
cls_grow(obj_array_grow1, obj_array_grow2, true)
ans(1).x

obj1.cls_insert([-1; -2; 0], [-3 -4 -2] ,1)
obj1.cls_insert([-1; -2; 0], [-3 -4 -2] ,1).cls_sort

%%

cls_grow(obj1, obj3, true)


obj_array5 = [obj1, obj2];
obj_array6 = [obj3, obj4];

a = cls_grow(obj_array6, obj_array5);
a(1)

%%

intersect([1,2,3,4], [2;3])

% test wrong interpolation
x_intr = [1, 2, 3, 4, 5, 6, 7, 5.5];
y_intr = x_intr - 2;
obj_intr = cls_line(x_intr, y_intr);
obj_intr.cls_interp([2;12])

a = [1,2;3,4]
sum(a, 1)
max(a)
repmat(ans, size(a,1), 1)

a = obj_intr.x(2:end)
class(a)
sort(obj_intr)

a = cls_line
a(1) = cls_line(x1, x1)
a(2) = cls_line(x1, x1)
a = reshape(a, [], 1)

