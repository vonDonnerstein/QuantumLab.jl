LIBINT_VERSION = "2.2.0-beta1"
lib_path = joinpath(Pkg.dir("QuantumLab"),"deps","usr","lib")

@linux_only begin
	cd(joinpath(Pkg.dir("QuantumLab"),"deps","usr","lib"))
	run(`ln -s libint2jl-$LIBINT_VERSION.so libint2jl.so`)
	run(`ln -s libint2-$LIBINT_VERSION.so libint2-QuantumLab.so`)
end

@osx_only begin
	cd(joinpath(Pkg.dir("QuantumLab"),"deps","usr","lib"))
	run(`ln -s libint2jl-$LIBINT_VERSION.dylib libint2jl.so`)
	run(`ln -s libint2-$LIBINT_VERSION.dylib libint2-QuantumLab.so`)
end

