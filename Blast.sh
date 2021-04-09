#Iniciamso sesion en tepeu para correr BLAST 
ssh -Y cmourra@tepeu.lcg.unam.mx 

# Nos dirigimos a la carpeta donde correremos el BLAST
cd /export/storage/users/cmourra/Clustering

# Corremos el BLAST con los parametros indicados
blastp -query ABC.faa -subject ABC.faa -outfmt 7 -max_hsps 1 -use_sw_tback -evalue 100 -out ABC3.faa.ali
