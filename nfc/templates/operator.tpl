class ${name}: public Tpetra::Operator<double,int,int> {

public:

${name}(
    const std::shared_ptr<const nosh::mesh> & mesh
    ):
  ${members_init}
  ${name}_light(mesh, $args)
{
}

virtual
~${name}()
{
}

void
apply(
    const Tpetra::MultiVector<double,int,int> & x,
    Tpetra::MultiVector<double,int,int> & y,
    Teuchos::ETransp mode,
    double alpha,
    double beta
    ) const
{
  this->${name}_light.apply(x, y, mode, alpha, beta);
  return;
}

Teuchos::RCP<const Tpetra::Map<int,int>>
getDomainMap() const
{
  return this->${name}_light.getDomainMap();
}

Teuchos::RCP<const Tpetra::Map<int,int>>
getRangeMap() const
{
  return this->${name}_light.getRangeMap();
}

protected:
private:
  ${members}
} // class ${name}


class ${name}_light: public Tpetra::Operator<double,int,int> {

public:

${name}_light(
    const std::shared_ptr<const nosh::mesh> & mesh
    ${arguments_light}
    ) :
  ${members_init_light}
{
}

virtual
~${name}_light()
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
${apply}
return;
}

Teuchos::RCP<const Tpetra::Map<int,int>>
getDomainMap() const
{
  return mesh_->map;
}

Teuchos::RCP<const Tpetra::Map<int,int>>
getRangeMap() const
{
  return mesh_->map;
}

protected:
private:
  ${members_light}
} // class ${name}_light
