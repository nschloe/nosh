
class ${name}: public Tpetra::Operator<double,int,int> {

public:

${name}(
    const std::shared_ptr<const Nosh::Mesh> & mesh,
    const Tpetra::Vector<double,int,int> & x0
    ) :
  mesh_(mesh),
  x0_(x0)
{
}

~${name}()
{
}

void
apply(
    const Tpetra::MultiVector<double,int,int> &x,
    Tpetra::MultiVector<double,int,int> &y,
    Teuchos::ETransp mode,
    double alpha,
    double beta
    ) const
{
TEUCHOS_TEST_FOR_EXCEPT_MSG(
    mode != Teuchos::NO_TRANS,
    "Only untransposed applies supported."
    );
TEUCHOS_TEST_FOR_EXCEPT_MSG(
    alpha != 1.0,
    "Only alpha==1.0 supported."
    );
TEUCHOS_TEST_FOR_EXCEPT_MSG(
    beta != 0.0,
    "Only beta==0.0 supported."
    )
${body}
return;
}

Teuchos::RCP<const Tpetra::Map<int,int>>
getDomainMap() const
{
  return x0_->getMap();
}

Teuchos::RCP<const Tpetra::Map<int,int>>
getRangeMap() const
{
  return x0_->getMap();
}

void
rebuild(
    const std::map<std::string, double> & params,
    const Tpetra::Vector<double,int,int> & x0
    )
{
  params_ = params;
  x0_ = x0;
  return;
}

protected:
private:
  const std::shared_ptr<const Nosh::Mesh> mesh_;
  Tpetra::Vector<double,int,int> x0_;
  std::map<std::string, double> params_;
} // class ${name}
