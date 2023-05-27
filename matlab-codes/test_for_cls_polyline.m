x1 = [1, 2, 3, 4, 5, 6, 7];
y1 = x1 - 2;

x2 = [2;3;4;5;6;7];
y2 = x2 * 2;

obj1 = cls_polyline(x1, y1);
obj2 = cls_polyline(x2, y2);
obj = [obj1, obj2];
numel(obj)
size(obj)

objarray = [obj1, obj2; obj2, obj1]

a = [-1,2,-3,4;5,-6,7,-8];
a(7)
b = NaN(1, numel(a))
for i = 1:numel(a)
  b(i) = a(i)
end
b = reshape(b, size(a))
a
sort(a)
[~, i] = sort(a)
ans = obj.x
unique([objarray.x])


obj1.x <= obj2.x




