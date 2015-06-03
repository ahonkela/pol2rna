#! /bin/zsh

MYDIR="$HOME/local/projects/synergy_rnaseq"
BITSEQ="$HOME/github/BitSeq/bitseq"

for k in `ls -1 ${MYDIR}/*.counts` ; do
    GENEFILE=${k/.counts/.genecounts}
    MEANFILE=${k/.counts/.genemeans}
    BASEFILE=`basename $k`
    echo Processing file $k
    "$BITSEQ" getGeneExpression -o "${GENEFILE}" -t "${MYDIR}/Homo_sapiens.GRCh37.68.cdna.new_ref.fixed.tr" -g "${MYDIR}/Homo_sapiens.GRCh37.68.cdna.new_ref.${BASEFILE}.genes" "$k"
    "$BITSEQ" getVariance -o "${MEANFILE}" "${GENEFILE}"
done
