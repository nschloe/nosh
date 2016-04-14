#ifndef NOSH_MODEL_H
#define NOSH_MODEL_H

// includes
#include <map>
#include <string>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_ParameterList.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_ModelEvaluatorDefaultBase.hpp>

namespace nosh
{
class model : public Thyra::ModelEvaluatorDefaultBase<double>
{
public:
  // f, jac, dfdp,
  // {
  //   {"linear solver package", "Belos"},
  //   {"method", "Pseudo Block CG"},
  //   // {"preconditioner", problem.prec}
  // }
  model (
      const std::shared_ptr<nosh::mesh> & mesh,
      const std::shared_ptr<nosh::fvm_operator> & f,
      const std::shared_ptr<nosh::fvm_operator> & jac,
      const std::shared_ptr<nosh::fvm_operator> & dfdp
    ):
    mesh_(mesh),
    f_(f),
    jac_(jac),
    dfdp_(dfdp),
    space_(Thyra::createVectorSpace<double>(Teuchos::rcp(mesh_->map())))
  {
    // Initialize the parameters
    const auto f_params = f_->get_scalar_parameters();
    const auto jac_params = jac_->get_scalar_parameters();
    const auto dfdp_params = dfdp_->get_scalar_parameters();

    std::map<std::string, double> all_params;
    all_params.insert(f_params.begin(), f_params.end());
    all_params.insert(jac_params.begin(), jac_params.end());
    all_params.insert(dfdp_params.begin(), dfdp_params.end());

    p_map_ = Teuchos::rcp(new Tpetra::Map<int,int>(
          all_params.size(),
          0,
          Teuchos::rcp(mesh_->comm)
          ));

    auto p_init = Thyra::createMember(this->get_p_space(0));
    p_names_ = Teuchos::rcp(new Teuchos::Array<std::string>(all_params.size()));
    int k = 0;
    for (auto it = all_params.begin(); it != all_params.end(); ++it) {
      (*p_names_)[k] = it->first;
      Thyra::set_ele(k, it->second, p_init());
      k++;
    }

    // set nominal values
    const Teuchos::RCP<const Tpetra::Vector<double,int,int>> initial_x =
      Teuchos::rcp(
        new Tpetra::Vector<double,int,int>(Teuchos::rcp(mesh_->map()))
        );
    const auto xxx = Thyra::createConstVector(initial_x, space_);
    nominal_values_.set_p(0, p_init);
    nominal_values_.set_x(xxx);
  }

  virtual
  ~model() {};

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_x_space() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(!space_.is_null());
#endif
    return space_;
  }

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_f_space() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(!space_.is_null());
#endif
    return space_;
  }

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_p_space(int l) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        l != 0,
        "Nosh can only deal with one parameter vector."
        );
    return Thyra::createVectorSpace<double>(p_map_);
  }

  virtual
  Teuchos::RCP<const Teuchos::Array<std::string> >
  get_p_names(int l) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        l != 0,
        "Nosh can only deal with one parameter vector."
        );
    return p_names_;
  }

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_g_space(int l) const
  {
    (void) l;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "Not implemented."
        );
    return Teuchos::null;
  }

  virtual
  Teuchos::ArrayView<const std::string>
  get_g_names(int j) const
  {
    (void) j;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "Not implemented."
        );
    return Teuchos::null;
  }

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  getNominalValues() const
  {
    return nominal_values_;
  }

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  getLowerBounds() const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "Not implemented."
        );
    return this->createInArgs();
  }

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  getUpperBounds() const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "Not implemented."
        );
    return this->createInArgs();
  }

  //virtual
  //Teuchos::RCP<Thyra::LinearOpWithSolveBase<double>>
  //create_W() const;

  virtual
  Teuchos::RCP<Thyra::LinearOpBase<double>>
  create_W_op() const
  {
    Teuchos::RCP<Tpetra::Operator<double,int,int>> op = Teuchos::rcp(jac_);
    return Thyra::createLinearOp(op, space_, space_);
  }

  virtual
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
  get_W_factory() const
  {
    // TODO make configurable
    Stratimikos::DefaultLinearSolverBuilder builder;

    auto p = Teuchos::rcp(new Teuchos::ParameterList);
    p->set("Linear Solver Type", "Belos");
    auto & belosList =
      p->sublist("Linear Solver Types")
      .sublist("Belos");
    //belosList.set("Solver Type", "MINRES");
    //belosList.set("Solver Type", "Pseudo Block GMRES");
    belosList.set("Solver Type", "Pseudo Block CG");

    auto & solverList =
      belosList.sublist("Solver Types")
      .sublist("Pseudo Block CG");
    solverList.set("Output Frequency", 1);
    solverList.set("Output Style", 1);
    solverList.set("Verbosity", 33);

    p->set("Preconditioner Type", "None");
    builder.setParameterList(p);

    auto lowsFactory = builder.createLinearSolveStrategy("");

    lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

    return lowsFactory;
  }

  virtual
  Teuchos::RCP<Thyra::PreconditionerBase<double>>
  create_W_prec() const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "Not implemented."
        );
    return Teuchos::null;
  }

  virtual
  void
  reportFinalPoint(
      const Thyra::ModelEvaluatorBase::InArgs<double> &finalPoint,
      const bool wasSolved
      )
  {
    (void) finalPoint;
    (void) wasSolved;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "Not implemented."
        );
  }

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  createInArgs() const
  {
    Thyra::ModelEvaluatorBase::InArgsSetup<double> in_args;

    in_args.setModelEvalDescription("Nosh problem");

    // TODO
    // We have *one* parameter vector with numParams_ parameters in it.
    in_args.set_Np(1);

    in_args.setSupports(IN_ARG_x, true);

    // for shifted matrix
    // TODO add support for operator shift
    in_args.setSupports(IN_ARG_alpha, true);
    in_args.setSupports(IN_ARG_beta, true);

    return in_args;
  }

protected:

  virtual
  Thyra::ModelEvaluatorBase::OutArgs<double>
  createOutArgsImpl() const
  {
    // TODO this entire method
    Thyra::ModelEvaluatorBase::OutArgsSetup<double> out_args;

    out_args.setModelEvalDescription("Nosh problem");

    out_args.set_Np_Ng(1, 0); // one parameter vector, no objective function

    out_args.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);

    // support derivatives with respect to all parameters;
    out_args.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,
        0,
        DerivativeSupport(DERIV_MV_BY_COL)
        );

    out_args.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
    // TODO
    // out_args.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec);

    return out_args;
  }

  virtual
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &in_args,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &out_args
      ) const
  {
    const double alpha = in_args.get_alpha();
    double beta = in_args.get_beta();

    // From packages/piro/test/MockModelEval_A.cpp
    if (alpha == 0.0 && beta == 0.0) {
      beta = 1.0;
    }
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(alpha, 0.0);
    TEUCHOS_ASSERT_EQUALITY(beta,  1.0);
#endif

    const auto & x_in = in_args.get_x();
#ifndef NDEBUG
    TEUCHOS_ASSERT(!x_in.is_null());
#endif
    // create corresponding tpetra vector
    auto x_in_tpetra =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getConstTpetraVector(
          x_in
          );

    // Dissect in_args.get_p(0) into parameter sublists.
    const auto & p_in = in_args.get_p(0);
#ifndef NDEBUG
    TEUCHOS_ASSERT(!p_in.is_null());
#endif

#ifndef NDEBUG
    // Make sure the parameters aren't NaNs.
    for (int k = 0; k < p_in->space()->dim(); k++) {
      TEUCHOS_ASSERT(!std::isnan(Thyra::get_ele(*p_in, k)));
    }
#endif

    // Fill the parameters into a std::map.
    const auto param_names = this->get_p_names(0);
    std::map<std::string, double> params;
    for (int k = 0; k < p_in->space()->dim(); k++) {
      params[(*param_names)[k]] = Thyra::get_ele(*p_in, k);
    }

    // compute F
    const auto & f_out = out_args.get_f();
    if (!f_out.is_null()) {

      auto f_out_tpetra =
        Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraVector(
            f_out
            );
      this->f_->set_parameters(params, {});
      this->f_->apply(
          *x_in_tpetra,
          *f_out_tpetra
          );
    }

    // Compute df/dp.
    const auto & derivMv = out_args.get_DfDp(0).getDerivativeMultiVector();
    const auto & dfdp_out = derivMv.getMultiVector();
    if (!dfdp_out.is_null()) {
      auto dfdp_out_tpetra =
        Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraMultiVector(
            dfdp_out
            );

      const int numAllParams = this->get_p_space(0)->dim();
      TEUCHOS_ASSERT_EQUALITY(
          numAllParams,
          dfdp_out_tpetra->getNumVectors()
          );
      // Compute all derivatives.
      this->dfdp_->set_parameters(params, {});
      for (int k = 0; k < numAllParams; k++) {
        this->dfdp_->apply(
            *x_in_tpetra,
            *dfdp_out_tpetra->getVectorNonConst(k)
            );
      }
    }

    // Fill Jacobian.
    const auto & W_out = out_args.get_W_op();
    if(!W_out.is_null()) {
      auto W_outT =
        Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraOperator(
            W_out
            );
      const auto & jac =
        Teuchos::rcp_dynamic_cast<nosh::fvm_operator>(W_outT, true);
      std::shared_ptr<const Tpetra::Vector<double,int,int>> x_std =
        Teuchos::get_shared_ptr(x_in_tpetra);
      const std::map<std::string, std::shared_ptr<const Tpetra::Vector<double, int, int>>> mp = {{"u0", x_std}};
      jac->set_parameters(params, mp);
    }

//     // Fill preconditioner.
//     const auto & WPrec_out = out_args.get_W_prec();
//     if(!WPrec_out.is_null()) {
//       auto WPrec_outT =
//         Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraOperator(
//             WPrec_out->getNonconstUnspecifiedPrecOp()
//             );
//       const auto & keoPrec =
//         Teuchos::rcp_dynamic_cast<nosh::keo_regularized>(WPrec_outT, true);
//       keoPrec->rebuild(
//           params,
//           *x_in_tpetra
//           );
//     }
    return;
  }

private:
  const std::shared_ptr<nosh::mesh> mesh_;
  const std::shared_ptr<nosh::fvm_operator> f_;
  const std::shared_ptr<nosh::fvm_operator> jac_;
  const std::shared_ptr<nosh::fvm_operator> dfdp_;

  Teuchos::RCP<const Tpetra::Map<int,int>> p_map_;
  Teuchos::RCP<Teuchos::Array<std::string> > p_names_;
  Thyra::ModelEvaluatorBase::InArgs<double> nominal_values_;

  const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> space_;
};
} // namespace nosh

#endif // NOSH_MODEL_H
