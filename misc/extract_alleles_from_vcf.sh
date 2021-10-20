doawk() {
	FILE=$1
	BASE=$2
	awk -F"\t" '!/#/ { print $1, ($2-1),$2, $4, $5}' "$FILE" > "${BASE}_alleles.bed"
}
export -f doawk
parallel doawk {} {.} ::: *_anc.vcf
doawk() {
	FILE=$1
	awk '{print >> FILENAME"."$1; close($1)}' "$FILE"
}

export -f doawk
parallel doawk {} ::: *_anc_alleles.txt

