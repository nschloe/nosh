// radius of the circle
R = 1.25 * 3.4;

// characteristic length of an edge
lcar = 1.0e-1;

Point(1) = {0, 0, 0, lcar};
Point(2) = {R, 0, 0, lcar};
Point(3) = {-R, 0, 0, lcar};
Point(4) = {0, R, 0, lcar};
Point(5) = {0, -R, 0, lcar};
Point(6) = {R*Cos(Pi/12.), R*Sin(Pi/12.), 0, lcar};
Point(7) = {R*Cos(-Pi/12.), R*Sin(-Pi/12.), 0, lcar};
Point(8) = {0.5*R, 0, 0, lcar};
Circle(1) = {6, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 7};
Line(5) = {7, 8};
Line(6) = {8, 6};
Line Loop(5) = {1, 2, 3, 4, 5, 6};
Plane Surface(6) = {5};

// Name the circle. This makes sure that edges and points aren't explicitly
// exported, only the triangles; those are what we need.
Physical Surface("my circle") = 6;
