// Gmsh file for a cube of edge length 10 with several ball-shaped cavities.

// Define characteristic lengths of the edges to be used for point definition.
// Different characteristic lengths could be used for the corner points of
// the cube and the around the cavity.
lcar1 = 0.5;

// Define the corner points of the cube.
Point(1) = {5.0,5.0,5.0,lcar1};
Point(2) = {5.0,5.0,-5.0,lcar1};
Point(3) = {5.0,-5.0,5.0,lcar1};
Point(4) = {5.0,-5.0,-5.0,lcar1};
Point(5) = {-5.0,5.0,5.0,lcar1};
Point(6) = {-5.0,5.0,-5.0,lcar1};
Point(7) = {-5.0,-5.0,5.0,lcar1};
Point(8) = {-5.0,-5.0,-5.0,lcar1};

// Define cube lines.
Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {1,5};
Line(4) = {2,4};
Line(5) = {2,6};
Line(6) = {3,4};
Line(7) = {3,7};
Line(8) = {4,8};
Line(9) = {5,6};
Line(10) = {5,7};
Line(11) = {6,8};
Line(12) = {7,8};

// Define cube surfaces
Line Loop(22) = {1,4,-6,-2};    Plane Surface(31) = {22};
Line Loop(23) = {1,5,-9,-3};    Plane Surface(32) = {23};
Line Loop(24) = {2,7,-10,-3};   Plane Surface(33) = {24};
Line Loop(25) = {4,8,-11,-5};   Plane Surface(34) = {25};
Line Loop(26) = {6,8,-12,-7};   Plane Surface(35) = {26};
Line Loop(27) = {9,11,-12,-10}; Plane Surface(36) = {27};

// Instead of using included files, we now use a user-defined function
// in order to carve some holes in the cube.
Function CheeseHole
  // In the following commands we use the reserved variable name
  // `newp', which automatically selects a new point number. This
  // number is chosen as the highest current point number, plus
  // one. (Note that, analogously to `newp', the variables `newc',
  // `news', `newv' and `newreg' select the highest number amongst
  // currently defined curves, surfaces, volumes and `any entities
  // other than points', respectively.)
  p1 = newp; Point(p1) = {x,  y,  z,  lcar1} ;
  p2 = newp; Point(p2) = {x+r,y,  z,  lcar1} ;
  p3 = newp; Point(p3) = {x,  y+r,z,  lcar1} ;
  p4 = newp; Point(p4) = {x,  y,  z+r,lcar1} ;
  p5 = newp; Point(p5) = {x-r,y,  z,  lcar1} ;
  p6 = newp; Point(p6) = {x,  y-r,z,  lcar1} ;
  p7 = newp; Point(p7) = {x,  y,  z-r,lcar1} ;

  c1 = newreg; Circle(c1) = {p2,p1,p7};
  c2 = newreg; Circle(c2) = {p7,p1,p5};
  c3 = newreg; Circle(c3) = {p5,p1,p4};
  c4 = newreg; Circle(c4) = {p4,p1,p2};
  c5 = newreg; Circle(c5) = {p2,p1,p3};
  c6 = newreg; Circle(c6) = {p3,p1,p5};
  c7 = newreg; Circle(c7) = {p5,p1,p6};
  c8 = newreg; Circle(c8) = {p6,p1,p2};
  c9 = newreg; Circle(c9) = {p7,p1,p3};
  c10 = newreg; Circle(c10) = {p3,p1,p4};
  c11 = newreg; Circle(c11) = {p4,p1,p6};
  c12 = newreg; Circle(c12) = {p6,p1,p7};

  // We need non-plane surfaces to define the spherical holes. Here we
  // use ruled surfaces, which can have 3 or 4 sides:
  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1};
  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2};
  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3};
  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4};
  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5};
  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6};
  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7};
  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8};

  // We then store the surface loops identification numbers in a list
  // for later reference (we will need these to define the final
  // volume):
  theloops[t] = newreg;
  Surface Loop(theloops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};

  thehole = newreg;
  // Don't define a volume for the ball -- it's supposed to be a cavity.
  Volume(thehole) = theloops[t] ;

Return

// Make some cavities.
x = -1.46317 ; y = -2.40772 ; z = 1.4293; r = 0.457792;
t = 1;
Call CheeseHole;

x = -1.49481 ; y = 3.86381 ; z = 2.32999; r = 0.423834;
t = 2;
Call CheeseHole;

x = 3.51252 ; y = -3.17091 ; z = 0.126542; r = 0.574615;
t = 3;
Call CheeseHole;

x = 1.80873 ; y = -1.3192 ; z = 0.449734; r = 0.815995;
t = 4;
Call CheeseHole;

x = 2.81175 ; y = -3.05816 ; z = 1.60204; r = 0.750539;
t = 5;
Call CheeseHole;

x = -1.64069 ; y = 0.642466 ; z = 1.57775; r = 0.758124;
t = 6;
Call CheeseHole;

x = -4.02354 ; y = 1.49486 ; z = -2.06116; r = 0.840947;
t = 7;
Call CheeseHole;

x = -1.23371 ; y = -4.12543 ; z = -3.45951; r = 0.725355;
t = 8;
Call CheeseHole;

x = -2.19689 ; y = -3.36984 ; z = 2.87567; r = 1.1617;
t = 9;
Call CheeseHole;

x = -3.4803 ; y = 2.57227 ; z = 0.0892171; r = 0.853531;
t = 10;
Call CheeseHole;

x = 3.05549 ; y = 1.74856 ; z = 2.26179; r = 0.82265;
t = 11;
Call CheeseHole;

x = -3.55543 ; y = 2.60709 ; z = 2.43619; r = 0.711617;
t = 12;
Call CheeseHole;

x = -1.70232 ; y = -4.56487 ; z = -0.28901; r = 0.42475;
t = 13;
Call CheeseHole;

x = -1.71015 ; y = 3.76747 ; z = -0.493561; r = 0.718423;
t = 14;
Call CheeseHole;

x = 0.839863 ; y = -3.57311 ; z = 4.26341; r = 0.409265;
t = 15;
Call CheeseHole;

x = 2.98551 ; y = 2.01366 ; z = -1.09315; r = 0.609814;
t = 16;
Call CheeseHole;

x = -4.32456 ; y = 4.29347 ; z = -4.05296; r = 0.627758;
t = 17;
Call CheeseHole;

x = -1.20392 ; y = -3.64129 ; z = -1.91227; r = 0.864392;
t = 18;
Call CheeseHole;

x = 3.75785 ; y = 4.00852 ; z = 4.17353; r = 0.790652;
t = 19;
Call CheeseHole;

x = -0.00189355 ; y = 1.9328 ; z = -4.14051; r = 0.83671;
t = 20;
Call CheeseHole;

// We define a physical volume for each hole:
//  Physical Volume (t) = thehole ;
//  Printf("Hole %g (center = {%g,%g,%g}, radius = %g) has number %g!",
//	 t, x, y, z, r, thehole) ;

// Define outer surface of the cube.
theloops[0] = newreg;
Surface Loop(theloops[0]) = {31,32,33,34,35,36};

// The volume of the cube, without the 5 holes, is now defined by 6
// surface loops: the first surface loop defines the exterior surface;
// the surface loops other than the first one define holes.  (Again,
// to reference an array of variables, its identifier is followed by
// square brackets):
v = newreg;
Volume(v) = {theloops[]};

// We finally define a physical volume for the elements discretizing
// the cube, without the holes (whose elements were already tagged
// with numbers 1 to 5 in the `For' loop):
Physical Volume (newv) = v;
