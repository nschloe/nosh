#include "bratu.hpp"
#include <nosh.hpp>

using list = std::map<std::string, boost::any>;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  const auto f = std::make_shared<bratu::f>(mesh);
  const auto jac = std::make_shared<bratu::jacobian>(mesh);
  const auto dfdp = std::make_shared<bratu::dfdp>(mesh);

  // const auto problem = bratu::bratu(mesh);

  // Create a model evaluator.
  // Can be used with everything in Trilinos that accepts such a thing.
  const auto model = std::make_shared<nosh::model>(
      mesh, f, jac, dfdp
      // {
      //   {"linear solver package", "Belos"},
      //   {"method", "Pseudo Block CG"},
      //   // {"preconditioner", problem.prec}
      // }
  );

  // Check out
  // https://trilinos.org/docs/dev/packages/nox/doc/html/parameters.html
  // for a full parameter description.
  nosh::nonlinear_solve(
      model,
      {
        {"method", "Newton"},
        {"NOX Status Test", list{
          {"Test 1", list {
            {"Test Type", "MaxIters"},
            {"Maximum Iterations", 1}
          }}
        }},
        {"NOX", list{
          {"Tolerance", 1.0e2},
          {"Max steps", 1},
          {"Printing", list{
           {"Output Information", list{
             {"Details", true},
             {"Outer Iteration", true},
             {"Outer Iteration Status Test", true},
             {"Inner Iteration", true},
             {"Linear Solver Details", true},
             {"Parameters", true},
             {"Warning", true},
             {"Debug", true},
             {"Test Details", true},
             {"Error", true},
             {"Stepper Iteration", true},
             {"Stepper Details", true},
             {"Stepper Parameters", true}
           }}
          }}
        }},
      }
      );
    // <ParameterList name="NOX Status Test" >
    //     <Parameter name="Test Type" type="string" value="Combo"/>
    //     <Parameter name="Number of Tests" type="int" value="4"/>
    //     <Parameter name="Combo Type" type="string" value="OR"/>
    //     <ParameterList name="Test 0">
    //         <Parameter name="Test Type" type="string" value="Combo"/>
    //         <Parameter name="Number of Tests" type="int" value="3"/>
    //         <Parameter name="Combo Type" type="string" value="AND"/>
    //         <ParameterList name="Test 0">
    //             <Parameter name="Test Type" type="string" value="NormF"/>
    //             <Parameter name="Norm Type" type="string" value="Two Norm"/>
    //             <Parameter name="Tolerance" type="double" value="1.0e-8"/>
    //             <Parameter name="Scale Type" type="string" value="Unscaled"/>
    //         </ParameterList>
    //         <ParameterList name="Test 1">
    //             <Parameter name="Test Type" type="string" value="NormUpdate"/>
    //             <Parameter name="Norm Type" type="string" value="Two Norm"/>
    //             <Parameter name="Tolerance" type="double" value="1.0e-8"/>
    //             <Parameter name="Scale Type" type="string" value="Unscaled"/>
    //         </ParameterList>
    //         <ParameterList name="Test 2">
    //             <Parameter name="Test Type" type="string" value="NormWRMS"/>
    //             <Parameter name="BDF Multiplier" type="double" value="1.0"/>
    //             <Parameter name="Tolerance" type="double" value="1.0"/>
    //             <Parameter name="Alpha" type="double" value="1.0"/>
    //             <Parameter name="Beta" type="double" value="0.5"/>
    //             <Parameter name="Relative Tolerance" type="double" value="1.0e-5"/>
    //             <Parameter name="Absolute Tolerance" type="double" value="1.0e-8"/>
    //         </ParameterList>
    //     </ParameterList>
    //     <ParameterList name="Test 1">
    //         <Parameter name="Test Type" type="string" value="MaxIters"/>
    //         <Parameter name="Maximum Iterations" type="int" value="25"/>
    //     </ParameterList>
    //     <ParameterList name="Test 2">
    //         <Parameter name="Test Type" type="string" value="NormF"/>
    //         <Parameter name="Norm Type" type="string" value="Two Norm"/>
    //         <Parameter name="Tolerance" type="double" value="1.0e-13"/>
    //         <Parameter name="Scale Type" type="string" value="Unscaled"/>
    //     </ParameterList>
    //     <ParameterList name="Test 3">
    //         <Parameter name="Test Type" type="string" value="FiniteValue"/>
    //     </ParameterList>
    // </ParameterList>

#if 0
  // nosh::nonlinear_solve(
  //     model, x,
  //     {
  //       {"method", "Newton"},
  //       {"Direction", list{
  //         {"method": "steepest descend"},
  //         {"forcing term alpha", 1.0}
  //       }},
  //       {"parameters", list{
  //         {"method", "polynomial"},
  //         {"max iters", 10},
  //         {"interpolation type", "quadratic"}
  //       }}
  //     }
  //     );

  /*
  nosh::parameter_continuation(
      model, x,
      {
        {"Stepper", list{
          {"Continuation Method", "Arc Length"},
          {"Contunuation Parameter", "mu"},
          {"Initial Value", 0.0},
          {"Max Value", 10.0},
          {"Max Nonlinear Iterations", 5},
        }};
        {"Step Size", list{
          {"Initial Step Size", 1.0e-2}
        }},
        {"jacobian solve", list{
          {"jacobian", jac},
          {"linear solver package", "Belos"},
          {"method", "Pseudo Block CG"},
          {"preconditioner", M}
        }},
        {"fpdp", dfdp}
      }
      );
      */

  nosh::write(x, "out.h5m");
#endif

  return EXIT_SUCCESS;
}
