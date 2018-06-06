using PyPlot
using PyCall
@pyimport cmocean.cm as cmaps

cstart = 97
cend = 122

all_char = [Char(i) for i = cstart:cend]
append!(all_char,' ')
Nchar = length(all_char)

function index_of_char(c)
    if c == ' '
        return Nchar
    else
        return Int(c)-cstart+1
    end
end


char_count = zeros(Int,Nchar)
pair_count = zeros(Int,Nchar,Nchar)

open("galaxy_wiki.txt") do file
    s = readstring(file)
    ns = length(s)

    for (c1,c2) in zip(s[1:end-1],s[2:end])
        c1 = lowercase(c1)    # convert to lower case if possible
        if c1 in all_char
            i1 = index_of_char(c1)
            char_count[i1] += 1

            if c2 in all_char
                i2 = index_of_char(c2)
                pair_count[i1,i2] += 1
            end
        end
    end
end

N = sum(char_count)
Npair = sum(pair_count)

p = char_count/N
jP = pair_count/Npair



## plotting

fig,ax = subplots()

x = 0:Nchar-1
matshow((jP'./p)',cmap=cmaps.thermal)

xticks(x,string.(all_char),fontsize=8)
yticks(x,string.(all_char),fontsize=8)

colorbar()

xlabel("P(x|y)")

savefig("cond_prob2.pdf")
close(fig)
