DefineConstant[
dx = {0.25, Name "Half widht of inner rectangle"}
];

size = 0.25;

Point(1) = {0, 0, 0, size};
Point(2) = {-1, 0, 0, size};
Point(3) = {0, -1, 0, size};
Point(4) = {1, 0, 0, size};
Point(5) = {0, 1, 0, size};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(1) = {1, 2, 3, 4};

Point(6) = {-dx, 0, 0, size/2};
Point(7) = {0, -dx, 0, size/2};
Point(8) = {dx, 0, 0, size/2};
Point(9) = {0, dx, 0, size/2};

Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
Line Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};
Physical Surface(1) = {1};

Physical Line(1) = {5, 6, 7, 8};
Physical Line(2) = {1, 2, 3, 4};