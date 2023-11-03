cd class_sphere
gfortran -c ../constants/dmdconstants.f90
gfortran -c ../class_event/class_event.f90
gfortran -c ./class_sphere.f90
gfortran -c ./test_sphere_class.f90
gfortran dmdconstants.o class_event.o class_sphere.o test_sphere_class.o -o testSphereClass.exc
./testSphereClass.exc
cd ..
