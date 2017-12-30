LIBINT_VERSION = "2.2.0-beta1"
lib_path = joinpath(dirname(@__FILE__),"usr","lib")


@static if is_linux() 
	try begin 
			run(`grep --ignore-case --silent 'avx' /proc/cpuinfo`)
			SIMD = "avx"
			cd(joinpath(dirname(@__FILE__),"usr","lib"))
			run(`ln -sf libint2jl-$LIBINT_VERSION-$SIMD.so libint2jl.so`)
			run(`ln -sf libint2-$LIBINT_VERSION-$SIMD.so libint2-QuantumLab.so`)
		end
	catch
		warn("Your system (CPU) does not support AVX.\n
							The Libint2 Library for Linux is only provided as AVX enabled version.\n
							Therefore a Libint2 stub is used.\n
							You should be able to start QuantumLab and perform calculation but
							with a drastically decreased performance.\n
							You can try to build your own version of libint2.so and libint2jl.so.
							As a reference take a look at the Makefile in QuantumLab/deps/usr/src.")
	end
end

@static if is_apple()
	try begin
			run(pipeline(`sysctl machdep.cpu.features`,`grep --ignore-case  --silent 'avx\|sse'`))
			SIMD = "avx"
			try run(pipeline(`sysctl machdep.cpu.features`,`grep --ignore-case --silent 'avx'`))
			catch 
				SIMD = "sse"
			end
			cd(joinpath(dirname(@__FILE__),"usr","lib"))
			run(`ln -sf libint2jl-$LIBINT_VERSION-$SIMD.dylib libint2jl.so`)
			run(`ln -sf libint2-$LIBINT_VERSION-$SIMD.dylib libint2-QuantumLab.so`)
		end
	catch
		warn("Your CPU does support neither SSE nor AVX
							The Libint2 Library for Linux is only provided as AVX enabled version.\n
							Therefore a Libint2 stub is used.\n
							You should be able to start QuantumLab and perform calculation but
							with a drastically decreased performance.\n
							You can try to build your own version of libint2.dylib and libint2jl.dylib.
							As a reference take a look at the Makefile in QuantumLab/deps/usr/src.")
	end
end

