#!/usr/bin/env julia
doc="
Parse jelo output into a hdf5 record. Call with no <infile> in order to read from stdin.

Usage:
  ./txt2h5.jl [options] <outfile>
  ./txt2h5.jl [options] <infile> <outfile>

Options:
  --help -h            Print this help.
"
using parse;
using DocOpt;
using HDF5
opts = docopt(doc);
if nothing == opts["<infile>"]
    t,d,c = texttoarrays([eachline(STDIN)...])
else
    t,d,c = filetoarrays(opts["<infile>"]);
end

h5open(opts["<outfile>"], "w") do outfile
    write(outfile, "t", t)
    write(outfile, "comments", c)
    for pi in 1:size(d)[1]
        for (i,label) in zip(1:6,["x","y","z","bx","by","bz"])
            write(outfile, "$(pi)/$label", d[pi,:,i])
        end
    end
end

