import csv

# Abrimos los documentos necesarios
file = open("../Downloads/ABC3.faa.ali","r")
new_file = open("../Documents/alignment_scores.csv","w")
new_file2 = open("../Documents/alignment_scores2.csv","w")
writer = csv.writer(new_file)
writer2 = csv.writer(new_file2)

# Listas que contendran la información que necesitamos
alignments = []
bitscores = []
values = []

# Recorrer el archivo del BLAST para obtener la información
for line in file:
    
    # Obtenemos el alineamiento y el bitscore
    if line.startswith("#"):
        continue
    line = line.replace("\n","")
    line = line.split("\t")
    seq1 = line[0]
    seq2 = line[1]
    bitscore = float(line[11])
    align = "{}_{}".format(seq1,seq2)
    val = (align,bitscore)

        
    # Rellenamos las listas con la informacion obtenida
    alignments.append(align)
    bitscores.append(bitscore)
    values.append(val)
  
    
# Generar diccionario que asocie el bitscore al alineamiento
align_scores = dict.fromkeys(alignments)
max_bit = max(bitscores)

# Generar los csv que seran interpretados por R
csv_rowlist = ["Alingment","Bitscore","Similarity","Dissimilarity"]
writer.writerow(csv_rowlist)
csv_rowlist2 = ["Sequence_1","Sequence_2","Bitscore","Similarity","Dissimilarity"]
writer2.writerow(csv_rowlist2)

# ahora se buscara generar el archivo para el score de disimilitud
Alignment = []
Bitscore = []
Similarity = []
Dissimilarity = []


# recorrer el diccionario para obtener el escore de disimilitud a partir del bitscore
for key in align_scores:
    for align in values:
        if key == align[0]:
            bits = align[1]
            similarity = bits / max_bit
            dissimilarity = 1 - similarity
            
            Alignment.append(key)
            Bitscore.append(bits)
            Similarity.append(similarity)
            Dissimilarity.append(dissimilarity)
            
            
            scores = (bits,similarity,dissimilarity)
            align_scores[key] = scores
            


# Escribir los archivos csv
for i in range(len(Alignment)):
    writer.writerow([Alignment[i],Bitscore[i],Similarity[i],Dissimilarity[i]])
    seqs = Alignment[i].split("_")
    Seq_1 = seqs[0]
    Seq_2 = seqs[1]
    writer2.writerow([Seq_1,Seq_2,Bitscore[i],Similarity[i],Dissimilarity[i]])

    
# cerrar los archivos 
file.close()
new_file.close()
new_file2.close()
