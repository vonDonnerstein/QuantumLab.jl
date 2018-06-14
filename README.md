
# <img src="http://schurkus.com/wp-content/uploads/2015/10/QuantumLab.png" alt=Q width="55">uantumLab

`Linux/Mac`: [![Build Status](https://travis-ci.org/vonDonnerstein/QuantumLab.jl.svg)](https://travis-ci.org/vonDonnerstein/QuantumLab.jl) `Win`: [![Build Status](https://ci.appveyor.com/api/projects/status/github/vonDonnerstein/QuantumLab.jl?branch=master&svg=true)](https://ci.appveyor.com/project/vonDonnerstein/quantumlab) `Test Coverage`: [![Coverage Status](https://coveralls.io/repos/vonDonnerstein/QuantumLab.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/vonDonnerstein/QuantumLab.jl?branch=master)

A Quantum Chemistry Package based on the Julia language.

## What it is: A narrative introduction
Think about real world chemistry: Huge reactors are needed to produce massive amounts of material. However, in order to develop the reactions it would be a terrible overhead if the chemists would have to tweak with the complex reactors all the time. That's the point of laboratories. Places, where everything one needs to conduct the experiments can be quickly and easily taken out of the drawer and flexibly combined.

Classical quantum codes (Gaussian, QChem, Turbomole, ...) are like reactors - highly performant but inflexible. QuantumLab on the other hand provides the "lab theoretician" - the ladies and gentlemen designing the theories and algorithms of quantum theory - with a laboratory full of tools. By making everything accessible from a simple but flexible julia input one can test theoretical and algorithmic ideas quickly and easily.

The interactively tested and developed input can also be written to text files and later `include`d for reuse.  We call these simple text files "experiment protocols" and they typically end in '.jl'.  Elaboration on and perfectioning of these protocols allows the development of new tools which can be easily added back to the bench drawer (possibly after publishing in a journal) - just commit them in the src/ folder of QuantumLab and send a pull request. QuantumLab is designed with the scientific spirit of collaboration and making ones results accessible to the community in mind. That's why we made this code open-source under the MIT licence: To make it easy to enable the lab community with your protocols.

## Why Julia?
Julia walks like Python, runs like C.

Meaning: Julia is an emerging language that allows for ease-of-use like e.g. Python in contrast to Fortran and still resembles Fortran and C (the classical HPC languages) rather well in performance. Furthermore, as it runs on top of the Julia just-in-time (JIT) compiler it gets cross-plattform capabilities for free.

## Getting Started
One of the most essential theories in quantum chemistry is Hartree-Fock (HF). So let's get started by considering how to compute the HF energy with QuantumLab. First, we'll need to bring QuantumLab up.
```jl
julia> using QuantumLab
```
Almost all functions in QuantumLab follow a naming convention: verbNoun1Noun2... where the Nouns become more specific left-to-right. This allows for easy tab-completion and also allows to search for functions rather easily.
```jl
?> Hartree
```
gives us the name of a function that seems to be just what we are looking for: `evaluateHartreeFock`. The help shows that it's an alias for evaluateSCF and tells us how to use it: `evaluateSCF(basis, geometry, initialGuess, electronNumber)`. We are also informed about the choices for types of the arguments. Let's start by choosing the molecule for which to compute the ground-state energy, say, the water molecule. Handily enough, we find the corresponding geometry file in .xyz format in the test folder of QuantumLab, but you can choose any other geometry for which you have the geometry file at hand. How can we read the Geometry in? Again tab-completion and `?>` are our friends.
```jl
julia> h2o = Geometry("test/h2o.xyz")
```
Next, we need to choose the basis set. If we don't have a corresponding specification file at hand, QuantumLab will obtain one from basis set exchange (please make sure to comply with their terms of use if you use this feature).
```jl
julia> sto3g = BasisSet("STO-3G")
```
For the initial guess we can simply take the ZeroGuess and the "number of closed shell orbitals" of water is 5. So now we can
```jl
julia> evaluateSCF(sto3g, h2o, ZeroGuess, 5)
```

## Where to go from here?
Documentation has a tendency to get out of sync with the codebase quickly. Deprecated documentation is, however, not only useless, but can even be misleading and thereby harmful. Consequently, all documentation must be part of the continuous integration cycle and kept consistent if it exists in several places. To avoid these problems, this readme is kept to the bare minimum. The actual full documentation happens in the form of docstrings and can be accessed by julia's help functionality, e.g.
```jl
?> LaplaceModule
```
The long-term plan is to combine those docstrings into a full-blown documentation here with the help of the Documenter.jl package. Also, all documentation will at some point contain example cases that are automatically checked with jldoctests. Until then, please refer to test/runtests.jl for an inspiring, however chaotic, collection of usage examples. As the aim is to test for almost all lines of code during continuous integration, you should be able to find an example case for pretty much all functionality there. If this temporary solution frustrates you as much as me and you want to help with the long term solution, your pull request is more than welcome. 
