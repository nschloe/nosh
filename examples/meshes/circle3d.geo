// radius of the circle
R = 1.25 * 3.4;

// characteristic length of an edge
lcar = 1.0e-1;

cl1 = 1;
Point(1) = {0, 0, 0, lcar};
Point(2) = {R, 0, 0, lcar};
Point(3) = {-R, 0, 0, lcar};
Point(4) = {0, R, 0, lcar};
Point(5) = {0, -R, 0, lcar};
Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
Circle(3) = {3, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// produce a flat cookie
out[] = Extrude{0,0,0.5}{ Surface{6}; };
