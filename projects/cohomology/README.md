# Constructing Generators of Cohomology Classes on Surfaces

This folder contains the code developed for the bachelor thesis "Constructing Generators of Cohomology Classes on Surfaces" by Josua Rieder under the supervision of Prof. Dr. Ralf Hiptmair. The code is written in C++20 and is documented using Doxygen.

We have decided to develop the project as a header-only library because we make frequent use of C++ templates.

## File list

* [**cohomology.h**](cohomology.h): Contains most of the actual computational infrastructure developed in this project.
* [**verification.h**](verification.h): Contains the verification routine.
* [**utility.h**](utility.h): Contains some smaller bits that are not directly related to the cohomological computations.
* [**visualization.h**](visualization.h): Contains functions to export the data structures to OBJ files. These have been used in conjunction with ParaView to create some of the visualizations seen in the thesis.
* [**test/test.cc**](test/test.cc): The main unit test file.
* [**test/examples.h**](test/examples.h): Contains some mesh generators and example scenarios.
* [**example.cc**](example.cc): A small example that illustrates the API.
