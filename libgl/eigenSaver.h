#include <Teuchos_RCP.hpp>

#include <LOCA_SaveEigenData_AbstractStrategy.H>
#include <LOCA_Parameter_SublistParser.H>

#include <Teuchos_ParameterList.hpp>

class EigenSaver : public LOCA::SaveEigenData::AbstractStrategy
{

 public:

// Actually suggested interface:
//    EigenSaver(
//      const Teuchos::RCP<LOCA::GlobalData>& global_data,
//      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
//      const Teuchos::RCP<Teuchos::ParameterList>& eigenParams      );

   EigenSaver( const std::string fileName );

   virtual ~EigenSaver();

   virtual NOX::Abstract::Group::ReturnType
   save ( Teuchos::RCP<std::vector<double> >       &evals_r,
          Teuchos::RCP<std::vector<double> >       &evals_i,
          Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_r,
          Teuchos::RCP<NOX::Abstract::MultiVector> &evecs_i  );

  private:
    std::string fileName_;

};