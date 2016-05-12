using BinDeps

@BinDeps.setup

libint2jl = library_dependency("libint2jl", runtime=true, validate=(name,handle)->true)
libint2   = library_dependency("libint2-alpha", runtime=true)
#= libeigen3 = library_dependency("libeigen3", aliases=["libeigen3-dev","eigen3"]) =#
#= libgmp3   = library_dependency("libgmp3-dev", aliases=["libgmp-dev","libgmp10-dev","libgmp10","gmp"]) =#

# package managers
#= provides(AptGet, Dict("libint2-dev" => libint2, =#
#= 									"libeigen3-dev" => libeigen3, =#
#= 									"libgmp3-dev" => libgmp3)) =#
#= provides(Yum,    Dict("libint2" => libint2, =#
#= 									"eigen3" => libeigen3, =#
#= 									"gmp" => libgmp3)) =#
#= provides(Sources, Dict( URI("http://de.archive.ubuntu.com/ubuntu/pool/universe/libi/libint2/libint2-dev_2.1.0~beta2-2_amd64.deb") => libint2)) =#


prefix = joinpath(BinDeps.depsdir(libint2jl))
lib_dir = BinDeps.libdir(libint2jl)
download_dir = BinDeps.downloadsdir(libint2jl)
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

provides(SimpleBuild,
				(@build_steps begin
					CreateDirectory(lib_dir)
					FileRule(joinpath(prefix,"usr","lib","libint2jl.so"),
					@build_steps begin
						ChangeDirectory(download_dir)
						`cp libint2jl.so $lib_dir`
					end)
					#= (@build_steps begin =#
					#= 	ChangeDirectory(lib_dir) =#
					#= 	`touch sd.test` =#
					#= 	 #1= push!(Libdl.DL_LOAD_PATH,lib_dir) =1# =#
					#= end) =#
				end), libint2jl, os = :Linux)
provides(SimpleBuild,
				(@build_steps begin
					CreateDirectory(lib_dir)
					FileRule(joinpath(prefix,"usr","lib","libint2-alpha.so.2"),
					@build_steps begin
						ChangeDirectory(download_dir)
						`cp libint2-alpha.so.2 $lib_dir`
					end)
					#= (@build_steps begin =#
					#= 	ChangeDirectory(lib_dir) =#
					#= 	`touch sd.test212` =#
					#= 	#1= Libdl.dlopen("libint2-alpha.so.2") =1# =#
					#= end) =#
				end), libint2, os = :Linux)


@BinDeps.install Dict(:libint2jl => :libint2jl,
							:libint2 => :libint2)
							#= :libeigen3 => :libeigen3, =#
							#= :libgmp3 => :libgmp3 =#
