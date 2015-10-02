using QuantumLab
using Base.Test

# write your own tests here
@test 1 == 1

# test readxyz
@test readxyz("h2o.xyz").atoms[2].element.symbol == "O"
@test -0.752 < readxyz("h2o.xyz").atoms[3].position.x < -0.750

