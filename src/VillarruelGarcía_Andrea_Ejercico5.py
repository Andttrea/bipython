"""
EJERCICIO 5: La cadena proteica más larga de cualquiera de los ORFs de una secuencia de DNA dada.

Este script encuentra y analiza todos los marcos de lectura abiertos (ORFs) en una secuencia de DNA,
tanto en la cadena directa como en la complementaria reversa, y determina cuál produce la proteína
más larga.

"""

from Bio.Seq import Seq

secuencia_input = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
secuencia = Seq(secuencia_input)
secuencia_complementaria = secuencia.complement()
secuencia_reversa = secuencia.reverse_complement()
rna_compelementaria = secuencia_complementaria.transcribe()
rna_reversa = secuencia_reversa.transcribe()
codones_paro = ["UAA", "UAG", "UGA"]  # declaramos los codones de paro

def orfs(rna_secuencia, codones_paro):
    """
    Encuentra todos los marcos de lectura abiertos (ORFs) en una secuencia de RNA.
    
    Parámetros:
    -----------
    rna_secuencia : Bio.Seq.Seq
        Secuencia de RNA en la que buscar los ORFs
    codones_paro : list
        Lista de codones de parada (stop codons)
    
    Retorna:
    --------
    list
        Lista de strings con todos los ORFs encontrados en los 3 marcos de lectura
    """
    orfs =  []
    rna = str(rna_secuencia)  # convertimos el objeto Seq a string
    longitud = len(rna)

    for frame in range(3):
        posiciones_inicio = []
        for indice in range(frame, longitud - 2, 3):
            if rna[indice:indice + 3] == "AUG":
                posiciones_inicio.append(indice)

        for inicio in posiciones_inicio:
            codones = []
            for indice in range(inicio, longitud - 2, 3):
                codon_actual = rna[indice:indice + 3]
                codones.append(codon_actual)
                if codon_actual in codones_paro and indice != inicio:
                    orfs.append("".join(codones))
                    break
            else:
                orfs.append("".join(codones))

    return orfs

def traducir(orfs):
    """
    Traduce una lista de ORFs (secuencias de RNA) a sus correspondientes secuencias proteicas.
    
    Parámetros:
    -----------
    orfs : list
        Lista de ORFs como strings de RNA
    
    Retorna:
    --------
    list
        Lista de secuencias proteicas traducidas (sin asteriscos de parada)
    """
    proteinas = []
    for orf in orfs:
        proteina = str(Seq(orf).translate()).rstrip("*")
        if proteina:
            proteinas.append(proteina)
    return proteinas

def proteina_mas_larga(proteinas):
    """
    Encuentra la proteína más larga de una lista de proteínas.
    
    Parámetros:
    -----------
    proteinas : list
        Lista de secuencias proteicas como strings
    
    Retorna:
    --------
    tuple
        Tupla con (proteína_más_larga, longitud_máxima)
    """
    seleccionada = ""
    longitud_maxima = 0
    for proteina in proteinas:
        longitud_actual = len(proteina)
        if longitud_actual > longitud_maxima:
            longitud_maxima = longitud_actual
            seleccionada = proteina
    return seleccionada, longitud_maxima

orfs_directos = orfs(rna_compelementaria, codones_paro)
orfs_reversos = orfs(rna_reversa, codones_paro)
todos_orfs = orfs_directos + orfs_reversos
proteinas = traducir(todos_orfs)
proteina_larga, longitud = proteina_mas_larga(proteinas)

print(f"Total de ORFs encontrados en la secuencia: {len(todos_orfs)}")
print(f"Proteina mas larga: {proteina_larga}, Longitud: {longitud} aminoacidos")
