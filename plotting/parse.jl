#!/usr/bin/env julia
module parse

export texttoarrays,filetoarrays

#we lisping now!
parseline(s) = map(float,
                   split(strip(s),
                         r"[ ,\,:,;]"))

function texttoarrays(lines)
    if length(lines) > 0
        #ugh.
        commentline(s)=ismatch(r" *#",s)
        commentlines = filter(line -> commentline(line), lines);
        #remove comented lines
        filter!(line ->!commentline(line),  lines);
        # I give up, this isn't python. I'll try to see
        # how I can do this as a Julian.
        t = zeros(length(lines));
        first = parseline(lines[1]);
        t[1] = first[1];
        #                   t             particle num  quantity
        d = zeros(Float64,length(lines), length(first)-1);
        d[1,:] = first[2:end];
        #fucking julia's type inference is stupid, or I am stupid.
        for (i,v) in enumerate(lines[2:end])
            v = parseline(v);
            t[i+1] = v[1];
            d[i+1,:] = v[2:end];
        end
        #shape is then (particle, time, quantity).
        t,reshape(d, div(length(first)-1,6),length(lines),length(first)-1,),commentlines
    else
        [],[[]],[];
    end
end
filetoarrays(file::IOStream)=texttoarrays(readlines(file))
filetoarrays(filename::String)=open(filetoarrays, filename)
end

