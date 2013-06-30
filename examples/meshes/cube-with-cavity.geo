// Gmsh file for a cube of edge length 10 with a ball-shaped cavity in the
// the middle of radius 1. The object is centered around (0,0,0) and aligned
// with the axes.

// Define characteristic lengths of the edges to be used for point definition.
// Different characteristic lengths could be used for the corner points of
// the cube and the around the cavity.
lcar1 = 0.3;

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
Function Cavity
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
x = 0.0; y = 0.0; z = 0.0; r = 1.0;
t = 1;
Call Cavity;

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
