#primer ejercicio
import csv
from Bio import Entrez
from Bio import SeqIO
from collections import Counter
email = "ariel.zambrano@est.ikiam.edu.ec"
archivo_entrada = "data/gstm.txt"
archivo_salida = "results/source.csv"

def source(email, archivo_entrada, archivo_salida):
    Entrez.email = email

    numeros_acceso = []
    with open(archivo_entrada, 'r') as archivo:
        for i in range(100):
            linea = archivo.readline().strip()
            numeros_acceso.append(linea)

    nombres_especies = []
    for numero_acceso in numeros_acceso:
        handle = Entrez.efetch(db="nucleotide", id=numero_acceso, rettype="gb", retmode="text")
        registro = SeqIO.read(handle, 'gb')
        partes_descripcion = registro.description.split()
        nombre_especie = " ".join(partes_descripcion[1:3])
        nombres_especies.append(nombre_especie)

    frecuencias_especies = Counter(nombres_especies)

    with open(archivo_salida, "w", newline="") as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow(["Especie", "Frecuencia"])
        for especie, frecuencia in frecuencias_especies.items():
            escritor.writerow([especie, frecuencia])

source(email, archivo_entrada, archivo_salida)

#------------------------------------------------------------------------------------------------------
import csv
from Bio import Entrez
from Bio import SeqIO
from collections import Counter
from Bio.SeqUtils import ProtParam
import matplotlib.pyplot as plt

def sequences():
    Entrez.email = 'ariel.zambrano@est.ikiam.edu.ec'
# Limitación de 100
    accession_numbers = []
    with open('data/gstm.txt', 'r') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            if i >= 100:
                break
            accession_numbers.append(row[0].strip())
# Reconocimiento de secuencia
    for accession_number in accession_numbers:
        handle = Entrez.efetch(db='nucleotide', id=accession_number, rettype='gb', retmode='text')
        record = SeqIO.read(handle, 'gb')
        dna_sequence = record.seq
# Traducción
        protein_sequence = dna_sequence.translate()
        peptides = protein_sequence.split('*')
# Selección de peptidos iniciadores con Metionina
        met_peptides = [peptide for peptide in peptides if peptide.startswith('M')]
# Peso molecular e inestabilidad
        peptide_data = []
        for peptide in met_peptides:
            analysis = ProtParam.ProteinAnalysis(str(peptide))
            molecular_weight = analysis.molecular_weight()
            instability_index = analysis.instability_index()
            peptide_data.append((peptide, molecular_weight, instability_index))
# Guardardado en CSV
            with open('results/glupeptides.csv', 'a') as f:
                f.write(f'{accession_number},{peptide},{molecular_weight},{instability_index}\n')

        # Crear el gráfico de dispersión
        mw_values = [data[1] for data in peptide_data]
        ii_values = [data[2] for data in peptide_data]

        plt.scatter (mw_values, ii_values, color='purple', edgecolors='black', alpha=0.8, s=50, marker='o')
        plt.xlabel('Peso Molecular')
        plt.ylabel('Índice de estabilidad')
        plt.title('Glucopeptidos')
        plt.savefig('results/glupeptidos.png')
        plt.close()

sequences()
