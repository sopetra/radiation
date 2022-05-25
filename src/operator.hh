#ifndef OPERATOR_HH
#define OPERATOR_HH
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

#include "bctype.hh"
//#include "exact.hh"
/** Lokalni operator za zadaću :
 *
 *      - div( grad u) = 0      u \Omega
 *                   u = g(x)   na \Gamma_D\subseteq\partial\Omega
 *        - grad u . n = j(x)   na \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * sa konformnim konačnim elementima svih tipova u svim dimenzijama
 *
 * \tparam BCType klasa koja indicira rubni uvjet
 */


template<typename BCType, typename FEM>
class DiffusionLOP : // derivacijska lista -- jakobijan i pattern računa PDELab
  public Dune::PDELab::NumericalJacobianApplyVolume  <DiffusionLOP<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianVolume       <DiffusionLOP<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<DiffusionLOP<BCType, FEM> >,
  public Dune::PDELab::NumericalJacobianBoundary     <DiffusionLOP<BCType, FEM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // Zastavice koje signaliziraju da na svakom elementu treba zvati:
  enum { doPatternVolume = true };  // metodu za računanje patterna (iz volumnih doprinosa)
  enum { doAlphaVolume = true };    // alpha_volume
  enum { doAlphaBoundary = true };  // alpha_boundary

  using  LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType ; // tip konačnog elementa ...-> tip lokalne baze


  DiffusionLOP(const BCType& bctype_, // boundary cond.type
               unsigned int intorder_=4) :
    bctype( bctype_ ), intorder( intorder_ ) // intorder = 4 kako bi program radio za k = 1,2 i 3
  {}

  // Računanje volumnog integrala
  // eg   = element (geometry): geometrija entiteta
  // lfsu = lokalni prostor funkcija za rješenje: lokalna verzija prostora konačnih elemenata, restrikcija na elem.
  // lfsv = lokalni prostor funkcija za test funkciju
  // x    = vektor koeficijenata rješenja
  // r    = lokalni rezidual: cilj je napuniti lokalni rezidual
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
      const int dim = EG::Geometry::coorddimension;
      using Gradient = Dune::FieldVector<double, dim>; // vektor gradijenata
      auto gt = eg.geometry().type();
      auto const & rule = Dune::QuadratureRules<double, dim>::rule(gt, intorder);
      for(auto qpoint : rule){
          auto xi = qpoint.position();
          auto const & phi = cache.evaluateFunction(xi, lfsu.finiteElement().localBasis()); // izračuna sve bazne fije u danoj integracijskoj točki
          double u = 0.0;
          //for(std::size_t i=0; i<lfsu.size(); ++i) // iteriramo po svim baznim funkcijama koje je cache izračunao
          //    u += x(lfsu, i) * phi[i]; // x nije običan vektor, dohvaća se iti stupanj slobode

          auto const & gradphihat = cache.evaluateJacobian(xi, lfsu.finiteElement().localBasis());//gradijent baznih fija na referentnom elementu, a nama treba na fizičkom -> množenje sa gradijentom geometrijskog preslikavanja
          auto const & jac = eg.geometry().jacobianInverseTransposed(xi);//transonirani inverz od jacobijana u kvadr. točki

          std::vector<Gradient> gradphi(lfsu.size());
          for(std::size_t i=0; i<lfsu.size(); ++i)
              jac.mv(gradphihat[i][0],gradphi[i]); // spremi se u gradphi; gradphi[i]=jac * gradphihat[i][0]

          Gradient gradu(0.0);
          for(std::size_t i=0; i<lfsu.size(); ++i)
              gradu.axpy(x(lfsu, i), gradphi[i]);//gradu += x(lfsu, i) * gradphi[i]; axpy napravi linearnu kombinaciju argumenata

          //auto glob_xi = eg.geometry().global(xi);
          //double f_coeff = RHS(glob_xi);
          //double a_coeff = react_coeff(glob_xi);
          double factor = qpoint.weight()*eg.geometry().integrationElement(xi);
          for(std::size_t i=0; i<lfsu.size(); ++i)
              r.accumulate(lfsu, i, factor*(gradu * gradphi[i]) );//  (gradu * gradphi[i]) = skalarni produkt 2 field vektora

      }
  }

  // integral po rubu
  // ig     = intersection (= stranica elementa)
  // lfsu_s = lokalni prostor funkcija na stranici za rješenje
  // lfsu_v = lokalni prostor funkcija na stranici za test funkciju
  // x_s    = vektor koeficijenata rješenja (na stranici)
  // r_s    = rezidual (na stranici)
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
      const int dim = IG::Geometry::coorddimension;
      auto gtface = ig.geometry().type();
      auto const & rule = Dune::QuadratureRules<double, dim-1>::rule(gtface, intorder);
      const double c = 1e-8, u_e = 25.00;

      for(auto qpoint : rule){
          auto xi = qpoint.position(); // xi je na referentnom elementu stranice

          if(bctype.isDirichlet(ig, xi))
              continue;

          auto local = ig.geometryInInside().global(xi); // preslikava xi sa referentnog elementa stranice na stranicu referentnog elementa elementa
          auto const & phi = cache.evaluateFunction(local, lfsu_s.finiteElement().localBasis()); // izračuna sve bazne fije u danoj integracijskoj točki

          double u = 0.0;
          for(std::size_t i=0; i<lfsu_s.size(); ++i) // iteriramo po svim baznim funkcijama koje je cache izračunao
              u += x_s(lfsu_s, i) * phi[i];

          double eps = 1E-6;
          auto glob_xi = ig.geometry().global(xi);
          double j_coeff;
          if(glob_xi[0]<eps || glob_xi[1]<eps || glob_xi[0]>1-eps || glob_xi[1]>1-eps) // izolirani
              j_coeff = 0.0;
          else
              j_coeff =  c * (std::pow(u_e + 273, 4) - std::pow(u + 273, 4)); // S-B zakon

          double factor = qpoint.weight()*ig.geometry().integrationElement(xi);

          for(std::size_t i=0; i<lfsu_s.size(); ++i)
              r_s.accumulate(lfsu_s, i, - j_coeff*phi[i]*factor );//  (gradu * gradphi[i]) = skalarni produkt 2 field vektora

      }


  }

private:
  BCType const & bctype; // identificira Dirichletovu granicu
  unsigned int intorder;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache; // jednom kad izračunamo bazne fije na lokalnom elementu, onda ne treba ponovo
                                                   // vrati vektor svih baznih fija u integracijskoj točki

};

#endif // OPERATOR_HH
