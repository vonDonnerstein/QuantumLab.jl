
# <img src="http://schurkus.com/wp-content/uploads/2015/10/QuantumLab.png" alt=Q width="55">uantumLab

[![Build Status](https://travis-ci.org/vonDonnerstein/QuantumLab.svg)](https://travis-ci.org/vonDonnerstein/QuantumLab)

A Quantum Chemistry Package based on the Julia language.

## What it is: A narrative introduction
Imagine real world chemistry: When one goes about producing massive amounts of material huge reactors are needed. But in order to develop the reactions it would be a terrible overhead if the chemists would have to tweak with the complex reactors all the time. That is where laboratories come in. Places, where everything one needs to conduct the experiments can be quickly and easily taken out of the drawer and flexibly combined.

While classical quantum codes (Gaussian, QChem, Turbomole, ...) are like reactors - highly performant but inflexible - QuantumLab provides the "lab theoretician" - the ladies and gentlemen designing the theories and algorithms of quantum theory - with a laboratory full of all the tools they might need. By making everything accessible from a simple but flexible julia input (which we term "experiment protocol") one can go about and test any theoretical and algorithmic idea very quickly. By elaborating on and perfectioning these protocols new tools are developed which can then easily be added back to the bench drawer. That's why we made this code open-source under the MIT licence: To make it easy to enable the lab community with your protocols.

## Why Julia?
Because it's an emerging language that allows for ease-of-use like e.g. Python in contrast to Fortran and still resembles Fortran and C (the classical HPC languages) rather well in performance. Furthermore, as it runs on top of the Julia just-in-time (JIT) compiler it gets cross-plattform capabilities for free.
