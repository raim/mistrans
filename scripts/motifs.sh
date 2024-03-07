

meme -minw 2 -nmotifs 3 seqcontext_fromto_Q\:G_long.fa -o QG
meme -minw 2 -nmotifs 3 seqcontext_fromto_T\:V_long.fa -o TV

momo modl -o QG_momo seqcontext_fromto_Q\:G_motif.fa 
momo modl -o TV_momo seqcontext_fromto_T\:V_motif.fa 
## NOTE: TM enriched; cf. https://www.nature.com/articles/srep40403
