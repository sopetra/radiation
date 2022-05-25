#ifndef DRIVER_HH
#define DRIVER_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>

#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "bctype.hh"
#include "operator.hh"

#include <memory>

// GV = GridView tip, gv je GridView.
template<class GV>
void driver(const GV& gv)
{
     // 1. Izbor konačnog elementa
    const int k = 2; // intorder = 4 u operator.hh
    using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV,double,double,k>;
    FEM fem(gv);

    // 2. Tip ograničenja
    using CONSTRAINTS = Dune::PDELab::ConformingDirichletConstraints;

    // 3. Tip vektora i matrica.
    using VBE = Dune::PDELab::ISTL::VectorBackend<>;
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MBE mbe(12);

    // 4. Prostor konačnih elemenata.
    using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CONSTRAINTS, VBE>;
    GFS gfs(gv, fem);

    // 5. Konstruiraj spremnika za ograničenja u prostoru konačnih elemenata.
    using CC = typename GFS::template ConstraintsContainer<double>::Type;
    CC cc;

    // 6. Asembliranje Dirichleteovih ograničenja.
    DirichletBdry bctype;
    Dune::PDELab::constraints(bctype, gfs, cc);

    // 7. Konstrukcija mrežnog operatora.
    using LOP = DiffusionLOP<DirichletBdry, FEM>;
    LOP lop(bctype);
    using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, double, double, double, CC, CC>;
    GO go(gfs, cc, gfs, cc, lop, mbe);

    // 8. Konstruirati vektor stupnjeva slobode
    using U = typename GO::Traits::Domain;
    U u(gfs, 0.0);

    // 9. Instancirati klasu koja daje Dirichletov rubni uvjet i ubaciti Dirichletove vrijednosti
    //    u vektor rješenja (stupnjeva slobode).
    BCExtension<GV> bce(gv);
    Dune::PDELab::interpolate(bce, gfs, u);
    //Dune::PDELab::set_nonconstrained_dofs(cc, 0.0, u);

    // 10. Izbor linearnog rješavača i prekondicionera.
    using LS = Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR;
    LS ls(5000);

    // 11. Instanciranje linearnog ili nelinearnog rješavača i rješavanje problema
    using SLP = Dune::PDELab::StationaryLinearProblemSolver<GO, LS, U>; //još jedan rješavač; za nelinearan zamijeniti sa Newton
    SLP slp(go, ls, u, 1E-10);
    slp.apply();

    //12. Grafički prikaz - kreiranje diskretnih mrežnih funkcija
    // Ako imamo egzaktno rješenje prikazujemo i njega i grešku.

    using DGF = Dune::PDELab::DiscreteGridFunction<GFS, U>;
    DGF u_dgf(gfs, u); //mrežna za u


    // 13. Grafički prikaz -- VTK writer
    Dune::VTKWriter<GV> writer(gv);
    //writer.addVertexData(u_dgf, "numeric"); // !!!Ne funkcionira jer u_dgf nije obični vektor
    writer.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter<DGF> >(u_dgf, "numeric") );
    writer.write("output");


}

#endif // DRIVER_HH
