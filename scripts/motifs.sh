
## TODO: link to specific installation

meme -minw 2 -nmotifs 3 seqcontext_fromto_Q\:G_long.fa -o meme/QG
meme -minw 2 -nmotifs 3 seqcontext_fromto_T\:V_long.fa -o meme/TV
meme -minw 2 -nmotifs 3 seqcontext_from_W_long.fa -o meme/W
meme -minw 2 -nmotifs 3 seqcontext_from_H_long.fa -o meme/H

momo modl -o momo/QG seqcontext_fromto_Q\:G_motif.fa 
momo modl -o momo/TV seqcontext_fromto_T\:V_motif.fa 
momo modl -o momo/W seqcontext_from_W_motif.fa 
momo modl -o momo/H seqcontext_from_H_motif.fa 
## NOTE: TM enriched; cf. https://www.nature.com/articles/srep40403

## TODO: clustalw on commanline

seq 
