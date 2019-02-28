DefineConstant[
dx = {0.25, Name "Half widht of inner rectangle"}
];

size = 0.25;

Point(1) = {0, 0, 0, size};
Point(2) = {1, 0, 0, size};
Point(3) = {1, 1, 0, size};
Point(4) = {0, 1, 0, size};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};


Point(5) = {0.5-dx, 0.5-dx, 0, size};
Point(6) = {0.5+dx, 0.5-dx, 0, size};
Point(7) = {0.5+dx, 0.5+dx, 0, size};
Point(8) = {0.5-dx, 0.5+dx, 0, size};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};
Physical Surface(1) = {1};

Physical Line(1) = {5, 6, 7, 8};
Physical Line(2) = {1, 2, 3, 4};