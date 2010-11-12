/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

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

#ifndef GINLA_STATUSTEST_FACTORY_H
#define GINLA_STATUSTEST_FACTORY_H

#include <Teuchos_RCP.hpp>
#include <LOCA_StatusTest_Abstract.H>
#include <LOCA_GlobalData.H>

#include "Ginla_StateTranslator_Virtual.h"

namespace Ginla {

namespace StatusTest {

class Factory
{

public:

  //! Constructor.
  Factory();

  //! Destructor.
  virtual ~Factory();

  //! Returns a status test set from a parameter list xml file.
  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildStatusTests( const std::string                                         & file_name,
                    const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & stateTranslator,
                    const Teuchos::RCP<const Teuchos::ParameterList>          & eigenInfo,
                    const Teuchos::RCP<const LOCA::GlobalData>                & globalData,
                    std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests = 0
                  ) const;

  //! Returns a status test set from a parameter list.
  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildStatusTests( Teuchos::ParameterList                                    & p,
                    const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & stateTranslator,
                    const Teuchos::RCP<const Teuchos::ParameterList>          & eigenInfo,
                    const Teuchos::RCP<const LOCA::GlobalData>                & globalData,
                    std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests=0
                  ) const;

private:

  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildComboTest( Teuchos::ParameterList                                    & p,
                  const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & stateTranslator,
                  const Teuchos::RCP<const Teuchos::ParameterList>          & eigenInfo,
                  const Teuchos::RCP<const LOCA::GlobalData>                & globalData,
                  std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests=0
                ) const;

  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildEnergyTest( Teuchos::ParameterList & p,
                   const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & stateTranslator
                 ) const;

  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildLoopTest( Teuchos::ParameterList& p,
                 const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & stateTranslator
               ) const;

  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildMaxAcceptedStepsTest( Teuchos::ParameterList & p
                           ) const;

  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildStabilityChangeTest( Teuchos::ParameterList & p,
                            const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo
                          ) const;

  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildTurnaroundTest( Teuchos::ParameterList & p
                     ) const;

  //! Checks if a tag is present in the param list and adds the test to the tagged_test std::map if true.  Returns true if a tag was present.
  bool checkAndTagTest( const Teuchos::ParameterList                   & p,
                        const Teuchos::RCP<LOCA::StatusTest::Abstract> & test,
                        std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests
                      ) const;
};

/*! \brief Nonmember helper function for the NOX::StatusTest::Factory.

\relates NOX::StatusTest::Factory

*/
Teuchos::RCP<LOCA::StatusTest::Abstract>
buildStatusTests( const std::string& file_name,
                  const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & stateTranslator,
                  const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo,
                  const Teuchos::RCP<const LOCA::GlobalData> & globalData,
                  std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests=0
                );

/*! \brief Nonmember helper function for the NOX::StatusTest::Factory.

\relates NOX::StatusTest::Factory

*/
Teuchos::RCP<LOCA::StatusTest::Abstract>
buildStatusTests( Teuchos::ParameterList& p,
                  const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & stateTranslator,
                  const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo,
                  const Teuchos::RCP<const LOCA::GlobalData> & globalData,
                  std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >* tagged_tests=0
                );


};

}

#endif // GINLA_STATUSTEST_FACTORY_H
