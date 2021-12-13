# exportar tot el path de les netcdf i intel:
export meu=/usr/local
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$meu/mpichIf/lib:$meu/netcdfIf/lib:$meu/netcdf-f-If/lib/:$meu/hdf5If/lib/:$meu/pnetcdfIf/lib

export pain=/opt/intel/oneapi/compiler/2021.4.0/linux
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$pain/lib:$pain/lib/x64:$pain/lib/emu:$pain/intel64_lin

export PATH=$PATH:$meu/mpichIf/bin:$meu/netcdfIf/bin:$meu/netcdf-f-If/bin:$pain/bin/intel64:$pain/bin 
