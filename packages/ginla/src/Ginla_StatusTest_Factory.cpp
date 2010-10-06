/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl|"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "Ginla_StatusTest_Factory.h"

// #include "Teuchos_TestForException.hpp"
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

// Concrete StatusTest Objects
#include <LOCA_StatusTest_Combo.H>

#include "Ginla_StatusTest_Energy.h"
#include "Ginla_StatusTest_Loop.h"
#include "Ginla_StatusTest_MaxAcceptedSteps.h"
#include "Ginla_StatusTest_ParameterLimits.h"
#include "Ginla_StatusTest_StabilityChange.h"
#include "Ginla_StatusTest_Turnaround.h"

// =============================================================================
Ginla::StatusTest::Factory::
Factory()
{}
// =============================================================================
Ginla::StatusTest::Factory::
~Factory()
{}
// =============================================================================
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::Factory::
buildStatusTests( const std::string& file_name ,
                  const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator,
                  const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo,
                  const Teuchos::RCP<const LOCA::GlobalData> & globalData,
                  std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests
                ) const
{   
  Teuchos::RCP<LOCA::StatusTest::Abstract> status_tests;
  
#ifdef HAVE_TEUCHOS_EXTENDED
  Teuchos::ParameterList param_list;
  Teuchos::updateParametersFromXmlFile("input.xml", &param_list);
  status_tests = this->buildStatusTests(param_list, stateTranslator, eigenInfo, globalData, tagged_tests);
#else
  std::string msg = "Error - Teuchos Extended Support must be enabled to use the xml reader for parameter lists.  Please rebuild the Trilinos Teuchos library with --enable-teuchos-extended in teh configure script.";
  TEST_FOR_EXCEPTION( true, std::logic_error, msg );
#endif

  return status_tests;
}
// =============================================================================
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::Factory::
buildStatusTests( Teuchos::ParameterList& p,
                  const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator,
                  const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo,
                  const Teuchos::RCP<const LOCA::GlobalData> & globalData,
                  std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests
                ) const
{ 
  Teuchos::RCP<LOCA::StatusTest::Abstract> status_test;

  std::string test_type = "???";

  if (Teuchos::isParameterType<std::string>(p, "Test Type"))
    test_type = Teuchos::get<std::string>(p, "Test Type");
  else {
    std::string msg = "Error - The \"Test Type\" is a required parameter in the Ginla::StatusTest::Factory!";
    TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }  

  if (test_type == "Combo")
    status_test = this->buildComboTest(p, stateTranslator, eigenInfo, globalData, tagged_tests);
  else if (test_type == "Energy")
    status_test = this->buildEnergyTest( p, stateTranslator );
  else if (test_type == "Loop")
    status_test = this->buildLoopTest( p, stateTranslator );
  else if (test_type == "Max accepted steps")
    status_test = this->buildMaxAcceptedStepsTest( p );
  else if (test_type == "Stability change")
    status_test = this->buildStabilityChangeTest( p, eigenInfo );
  else if (test_type == "Turnaround")
    status_test = this->buildTurnaroundTest( p );
  else {
    std::ostringstream msg;
    msg << "Error - the test type \"" << test_type << "\" is invalid!";
    TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  this->checkAndTagTest(p, status_test, tagged_tests);

  return status_test;
}
// =============================================================================
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::Factory::
buildComboTest( Teuchos::ParameterList& p,
                const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator,
                const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo,
                const Teuchos::RCP<const LOCA::GlobalData> & globalData,
                std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests
              ) const
{ 

  int number_of_tests = Teuchos::get<int>(p, "Number of Tests");
  
  std::string combo_type_string = Teuchos::get<std::string>(p, "Combo Type");
  LOCA::StatusTest::Combo::ComboType combo_type;
  if (combo_type_string == "AND")
    combo_type = LOCA::StatusTest::Combo::AND;
  else if (combo_type_string == "OR")
    combo_type = LOCA::StatusTest::Combo::OR;
  else
  {
    std::string msg = 
      "Error - The \"Combo Type\" must be \"AND\" or \"OR\"!";
    TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  
  Teuchos::RCP<LOCA::StatusTest::Combo> combo_test =
    Teuchos::rcp(new LOCA::StatusTest::Combo(combo_type, globalData));
  
  for (int i=0; i < number_of_tests; ++i)
  {
    ostringstream subtest_name;
    subtest_name << "Test " << i;
    Teuchos::ParameterList& subtest_list = p.sublist(subtest_name.str(), true);
    
    Teuchos::RCP<LOCA::StatusTest::Abstract> subtest =
      this->buildStatusTests(subtest_list, stateTranslator, eigenInfo, globalData, tagged_tests);
    
    combo_test->addStatusTest(subtest);
  }
  
  return combo_test;
}
// =============================================================================
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::Factory::
buildEnergyTest( Teuchos::ParameterList & p,
                 const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator
               ) const
{
  double maxFreeEnergy = Teuchos::get<double>( p, "Maximum free energy" );
  
  Teuchos::RCP<Ginla::StatusTest::Energy> status_test =
    Teuchos::rcp( new Ginla::StatusTest::Energy( stateTranslator, maxFreeEnergy ) );

  return status_test;
}
// =============================================================================
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::Factory::
buildLoopTest( Teuchos::ParameterList & p,
               const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator
             ) const
{  
  Teuchos::RCP<Ginla::StatusTest::Loop> status_test =
    Teuchos::rcp( new Ginla::StatusTest::Loop( stateTranslator ) );

  return status_test;
}
// =============================================================================
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::Factory::
buildMaxAcceptedStepsTest( Teuchos::ParameterList & p
                         ) const
{
  int maxAcceptedSteps = Teuchos::get<int>( p, "Max accepted steps" );
  
  Teuchos::RCP<Ginla::StatusTest::MaxAcceptedSteps> status_test =
    Teuchos::rcp( new Ginla::StatusTest::MaxAcceptedSteps( maxAcceptedSteps ) );

  return status_test;
}
// =============================================================================
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::Factory::
buildStabilityChangeTest( Teuchos::ParameterList & p,
                          const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo
                        ) const
{
  int stabilityChangeThreshold = Teuchos::get<int>( p, "Stability change threshold" );

  Teuchos::RCP<Ginla::StatusTest::StabilityChange> status_test =
    Teuchos::rcp( new Ginla::StatusTest::StabilityChange( eigenInfo,
                                                          stabilityChangeThreshold )
                );

  return status_test;
}
// =============================================================================
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::Factory::
buildTurnaroundTest( Teuchos::ParameterList & p
                   ) const
{
  Teuchos::RCP<Ginla::StatusTest::Turnaround> status_test =
    Teuchos::rcp( new Ginla::StatusTest::Turnaround() );

  return status_test;
}
// =============================================================================
bool
Ginla::StatusTest::Factory::
checkAndTagTest( const Teuchos::ParameterList& p,
                 const Teuchos::RCP<LOCA::StatusTest::Abstract>& test,
                 std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests
               ) const
{
  if ( (Teuchos::isParameterType<std::string>(p, "Tag")) && (tagged_tests != NULL) )
  {
    (*tagged_tests)[Teuchos::getParameter<std::string>(p, "Tag")] = test;
    return true;
  }

  return false;
}
// =============================================================================
// Nonmember function
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::
buildStatusTests( const std::string& file_name,
                  const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator,
                  const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo,
                  const Teuchos::RCP<const LOCA::GlobalData> & globalData,
                  std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests
                )
{
  Ginla::StatusTest::Factory factory;
  return factory.buildStatusTests(file_name, stateTranslator, eigenInfo, globalData, tagged_tests);
}
// =============================================================================
// Nonmember function
Teuchos::RCP<LOCA::StatusTest::Abstract>
Ginla::StatusTest::
buildStatusTests( Teuchos::ParameterList& p,
                  const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator,
                  const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo,
                  const Teuchos::RCP<const LOCA::GlobalData> & globalData,
                  std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests
                )
{
  Ginla::StatusTest::Factory factory;
  return factory.buildStatusTests(p, stateTranslator, eigenInfo, globalData, tagged_tests);
}
// =============================================================================
