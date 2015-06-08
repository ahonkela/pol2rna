function classify_seqname(name)
    if name[1:4] == "ENST"
        return 1
    elseif name[1:4] == "ENSG"
        return 2
    else
        println("Unknown sequence: $name")
        return 4
    end
end

times=["0000", "0005", "0010", "0020", "0040", "0080", "0160", "0320", "0640", "1280"]

datadir = "/cs/taatto/group/mlb/synergy/data/2012-03_RNA/Aligned_to_GRCh37_68/alignments"
samtools = homedir()*"/bin/samtools"

for myind in ARGS
    myi = int(myind)
    fp0 = open(`$samtools view $datadir/MCF7_t$(times[myi])min_RNA_v68_2012-03.bam` |> `cut -f 1,3` |> `grep -v 1:N:0:`, "r")
    fp = fp0[1]

    println("Processing file MCF7_t$(times[myi])min_RNA_v68_2012-03.bam")
    l = readline(fp)
    linesdone = 1
    t = split(l, '\t')
    myread = t[1]
    myflags = classify_seqname(t[2])
    counts = zeros(Int64, 7)
    for l in eachline(fp)
        linesdone += 1
        if linesdone % 1000000 == 0
            println("$linesdone lines done, counts: $counts")
        end
        t = split(l, '\t')
        if t[1] == myread  # new alignment for the same read
            myflags |= classify_seqname(t[2])
        else   # new read
            counts[myflags] += 1
            myread = t[1]
            myflags = classify_seqname(t[2])
        end
    end
    counts[myflags] += 1

    println("Final counts: $counts")
    println("mRNA, pre-mRNA, both")
    close(fp)
    outf = open("$datadir/MCF7_t$(times[myi])min_RNA_v68_2012-03.summarycounts", "w")
    println(outf, "mRNA pre-mRNA both error")
    println(outf, "$(counts[1]) $(counts[2]) $(counts[3]) $(sum(counts[4:7]))")
    close(outf)
end
