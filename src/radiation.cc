#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <dune/common/parallel/mpihelper.hh> 
#include<dune/common/parametertreeparser.hh>

#include <dune/grid/uggrid.hh>             // Koristimo UGGrid
#include <dune/grid/io/file/gmshreader.hh> // GmshReader klasa
#include <dune/grid/common/gridinfo.hh>

#include "driver.hh"

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);   
    const int dim = 3;
    using GridType = Dune::UGGrid<dim>;
    using GridView = typename GridType::LeafGridView;
    // Učitaj mrežu
    std::unique_ptr<GridType> gridptr = Dune::GmshReader<GridType>::read("src_dir/transfinite-3D-simplex.msh");

    int nref = 0;

    auto gv = gridptr->leafGridView();
    std::cout << " Učitavanje mreže je gotovo.\n";


    driver(gv);

    return 0;

}
