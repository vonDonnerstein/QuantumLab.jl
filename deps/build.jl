LIBINT_VERSION = "2.2.0-beta1"
lib_path = joinpath(Pkg.dir("QuantumLab"),"deps","usr","lib")


@linux_only begin
	try run(pipeline(`grep --ignore-case 'avx\|sse' /proc/cpuinfo`))
	catch ; error("Your CPU does support neither SSE nor AVX")
	end
	SIMD = "avx"
	try run(`grep 'avx' --ignore-case /proc/cpuinfo`)
	catch ; SIMD = "sse"
	end
	cd(joinpath(Pkg.dir("QuantumLab"),"deps","usr","lib"))
	run(`ln -s libint2jl-$LIBINT_VERSION-$SIMD.so libint2jl.so`)
	run(`ln -s libint2-$LIBINT_VERSION-$SIMD.so libint2-QuantumLab.so`)
end

@osx_only begin
	try run(pipeline(`sysctl machdep.cpu.features`,`grep --ignore-case 'avx\|sse'`))
	catch ; error("Your CPU does support neither SSE nor AVX")
	end
	SIMD = "avx"
	try run(pipeline(`sysctl machdep.cpu.features`,`grep --ignore-case 'avx'`))
	catch ; SIMD = "sse"
	end
	cd(joinpath(Pkg.dir("QuantumLab"),"deps","usr","lib"))
	run(`ln -s libint2jl-$LIBINT_VERSION-$SIMD.dylib libint2jl.so`)
	run(`ln -s libint2-$LIBINT_VERSION-$SIMD.dylib libint2-QuantumLab.so`)
end

