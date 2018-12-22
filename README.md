In my case(macOS High Sierra 10.13.6, gcc (MacPorts gcc5 5.5.0_1) 5.5.0), I installed gsl with brew. it will get installed in /usr/local/Celler/gsl.
Note: you would see below if you do compile with normal step
```
Undefined symbols for architecture x86_64:
  "_gsl_sf_bessel_J0", referenced from:
      _main in besel_exam-72d841.o
ld: symbol(s) not found for architecture x86_64
clang: error: linker command failed with exit code 1
```
Instead of the normal step, you can compile your C file firstly and link it to the library like below


```
brew install gsl
gcc -Wall -I/usr/local/Cellar/gsl/2.5/include -c schrodinger.c -o a.out
gcc -L/usr/local/Cellar/gsl/2.5/lib a.out -lgsl -o output.out
./output.out
```

Of course, /usr/local/Cellar/gsl/ should be replaced for the path where you installed gsl.
Warning: you would see error message such

```
dyld: lazy symbol binding failed: Symbol not found: _cblas_dnrm2
  Referenced from: /usr/local/opt/gsl/lib/libgsl.23.dylib
  Expected in: flat namespace

dyld: Symbol not found: _cblas_dnrm2
  Referenced from: /usr/local/opt/gsl/lib/libgsl.23.dylib
  Expected in: flat namespace

Abort trap: 6
```

if you got above error message, you could link another library to your application like this
```
gcc -L/usr/local/Cellar/gsl/2.5/lib a.out -lgslcblas -lgsl -o output.out
```

Finally you can enjoy GSL!