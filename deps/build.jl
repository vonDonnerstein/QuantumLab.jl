using BinDeps

@BinDeps.setup

libint2jl = library_dependency("libint2jl.so", validate=(name,handle)->true)
libint2   = library_dependency("libint2.so")

prefix = joinpath(BinDeps.depsdir(libint2jl))
lib_dir = BinDeps.libdir(libint2jl)
bin_dir = BinDeps.downloadsdir(libint2jl)
#= download_dir = joinpath(prefix,"downloads") =#
#= provides(SimpleBuild, =#
#= 				(@build_steps begin =#
#= 					CreateDirectory(download_dir) =#
#= 					(@build_steps begin =#
#= 						ChangeDirectory(download_dir) =#
#= 						FileRule(joinpath(prefix,"usr","lib","libint2.so"), @build_steps begin =#
#= 							`wget http://de.archive.ubuntu.com/ubuntu/pool/universe/libi/libint2/libint2-dev_2.1.0~beta2-2_amd64.deb` =#
#= 							`ar -vx libint2-dev_2.1.0~beta2-2_amd64.deb` =#
#= 							`tar -xvzf data.tar.gz -C $prefix` =#
#= 							`rm -f debian-binary control.tar.gz data.tar.xz` =#
#= 						end) =#
#= 					end) =#
#= 				end), libint2, os = :Linux) =#

#Sets symbolic link for precompiled libint2.so and libint2jl.so in deps/bin directory for Linux
provides(SimpleBuild,
				(@build_steps begin
					CreateDirectory(lib_dir)
					FileRule(joinpath(prefix,"usr","bin","libint2jl.so"),
					@build_steps begin
						ChangeDirectory(bin_dir)
						`ln -s libint2jl.so $bin_dir`
					end)
				end), libint2jl, os = :Linux)
provides(SimpleBuild,
				(@build_steps begin
					CreateDirectory(lib_dir)
					FileRule(joinpath(prefix,"usr","lib","libint2.so"),
					@build_steps begin
						ChangeDirectory(bin_dir)
						`ln -s libint2-2.2.0-beta1.so $(joinpath(lib_dir,"libint2.so"))`
					end)
				end), libint2, os = :Linux)

#Sets symbolic link for precompiled libint2.so and libint2jl.so in deps/bin directory for OSX (Darwin)
provides(SimpleBuild,
				(@build_steps begin
					CreateDirectory(lib_dir)
					FileRule(joinpath(prefix,"usr","lib","libint2jl.so"),
					@build_steps begin
						ChangeDirectory(bin_dir)
						`ln -s libint2jl.dylib $(joinpath(lib_dir,"libint2jl.so"))`
					end)
				end), libint2jl, os = :Darwin)
provides(SimpleBuild,
				(@build_steps begin
					CreateDirectory(lib_dir)
					FileRule(joinpath(prefix,"usr","lib","libint2.so"),
					@build_steps begin
						ChangeDirectory(bin_dir)
						`ln -s libint2-2.2.0-beta1.dylib $(joinpath(lib_dir,"libint2.so"))`
					end)
				end), libint2, os = :Darwin)


@BinDeps.install Dict(:libint2jl => :libint2jl,
							:libint2 => :libint2)
