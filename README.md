# Grobner basis calculation using FGB
This package provides a high-level ocaml limited interface to the [Arb library](https://arblib.org/index.html). For more information on functionality see the arb documentation.

## Dependencies
This package requires arb and flint as well as gmp to be installed for this package to build.

For ocaml opam is also helpful. To install these dependencies on Ubuntu run:

`sudo apt-get install opam libgmp-dev libflint-dev libflint-arb-dev`

To initialize opam run:

`opam init`

To install dune run:

`opam install dune`

### Optional dependency Zarith
This package contains two public libraries. arb, and optionally arb.zarith. In the arb library the only way to interface to the arb integer type is by the use of strings. The conversion from these strings to gmp values gives a slight overhead. arb.zarith has all the functionality as arb except it also provides conversions from zarith integers to arb integers directly. arb.zarith will be built and installed if you have zarith installed. Otherwise, only arb will be installed.

To install zarith run:

`opam install zarith`

## Building
To build the libraries run:

`dune build`

## Installing
To install the library as a dev version through opam run:

`opam install .`

## Usage
If the package is installed, documentation can be built with `odig odoc` and viewed with `odig doc`. Or documentation can be built locally with `dune build @doc`, and then viewed in _build/default/_doc/_html.

bin/main.ml gives some examples of using the library. To run use `dune exec bin/main.exe`