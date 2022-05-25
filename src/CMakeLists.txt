add_executable("radiation" radiation.cc driver.hh operator.hh bctype.hh)
target_link_dune_default_libraries("radiation")
