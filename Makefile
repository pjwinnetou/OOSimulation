all: run_macro

clean:
	rm run_macro

run_macro: run_macro.cc
	g++ run_macro.cc ../pythia8release/releases/lib/libpythia8.a -o run_macro -I../pythia8release/releases/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,../pythia8release/releases/lib -ldl `root-config --libs --cflags` -Wdelete-non-virtual-dtor

run_macro_MB: run_macro_MB.cc
	g++ run_macro_MB.cc ../pythia8release/releases/lib/libpythia8.a -o run_macro_MB -I../pythia8release/releases/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,../pythia8release/releases/lib -ldl `root-config --libs --cflags` -Wdelete-non-virtual-dtor
