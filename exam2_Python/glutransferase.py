from Bio import Entrez
from Bio import SeqIO
import csv

def source(especie):
    Entrez.email = "tu_correo_electronico@example.com"
    handle = Entrez.esearch(db="nucleotide", term=especie, retmax=10)
    record = Entrez.read(handle)
    ids = record["IdList"]
    
    especies_contadas = {}
    
    for id in ids:
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        organism = record.annotations["organism"]
        
        especies_contadas[organism] = especies_contadas.get(organism, 0) + 1
    
    # Guardar los resultados en un archivo CSV
    with open("results/source.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Especie", "Frecuencia"])
        
        for especie, frecuencia in especies_contadas.items():
            writer.writerow([especie, frecuencia])
    
    return especies_contadas
