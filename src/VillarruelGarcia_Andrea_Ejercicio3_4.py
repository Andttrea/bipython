#Usando el concepto de Herencia  una subclase de la clase gen, llamada tRNA, y otra clase llamada RNA no codificante. 
#  Luego deriva de tRNA otra subclase llamada proteina. 
#ejercicio 4: Genera una función longitud para la clase tRNA y una función longitud para la clase proteína.  
#La primera regresa el numero de nucleótidos, la segunda, el número de nucleótidos y de aminoácidos.


from random import choice

class Gen:
    nombre = ""
    inicio = 0
    final = 100
    secuencia = []

    def __init__(self, nombre, pos_inicio, pos_final, secuencia):
        self.nombre = nombre
        self.inicio = pos_inicio
        self.final = pos_final
        self.secuencia = secuencia

    def longitud(self):
        return self.final - self.inicio


    def gc_content(self):
        gc_count = sum(1 for base in self.secuencia if base in "GC")
        return gc_count / len(self.secuencia) if self.secuencia else 0

#ejercicio 3
class tRNA(Gen):
    tiene_anticodon = True 
    
    def __init__(self, nombre, pos_inicio, pos_final, secuencia, aminoacido="", anticodon=""): #overrinding 
        Gen.__init__(self, nombre, pos_inicio, pos_final, secuencia) #aqui llamamos al constructor de la clase padre 
        self.aminoacido = aminoacido
        self.anticodon = anticodon
        self.cargado = False  #Es para checar si el aminoacido esta cargado

    def cargar_aminoacido(self, aminoacido): 
        self.aminoacido = aminoacido
        self.cargado = True
        return f"El tRNA {self.nombre} esta cargado con el aminoacido {self.aminoacido}"
    
    def verificar_anticodon(self,codon):
        complementos = {"A":"U", "U":"A", "C":"G", "G":"C"}
        if len(codon) == 3 and len(self.anticodon) == 3:
            anticodon_complementario = "".join(complementos[base] for base in codon)
            return anticodon_complementario == self.anticodon
        return False
    
#ejercicio 4
    def longitud(self):
        numero_nucleoticos = len(self.secuencia)
        return f'El tRNA tiene {numero_nucleoticos} nucleotidos.'
    
class RNA_no_codificante(Gen):
    es_codificante = False  
        
    def __init__(self, nombre, pos_inicio, pos_final, secuencia, funcion_reguladora = "desconocida"):
        Gen.__init__(self, nombre, pos_inicio, pos_final, secuencia)
        self.funcion_reguladora = funcion_reguladora
    
    def formar_estructura_secundaria(self):
        estructuras_nocRNA = ["hairpin", "bulge", "internal loop"]
        return choice(estructuras_nocRNA)
        
        
class Proteina(tRNA):
    polipeptido = True
    
    def __init__(self, nombre, pos_inicio, pos_final, secuencia, aminoacido="", anticodon="", estructura_3d = "desconocida"):
        tRNA.__init__(self, nombre, pos_inicio, pos_final, secuencia, aminoacido, anticodon)
        self.plegada = False
        self.estructura_3d = estructura_3d

    def plegar_proteina(self, temperatura =37):
        if temperatura >= 20 and temperatura <= 40:
            self.plegada = True
            estructuras = ["alfa-helice", "beta plegada"]
            self.estructura_3d = choice(estructuras)
            return f"Proteina plegada en estructura {self.estructura_3d}"
        else:
            self.plegada = False
            self.estructura_3d = "desconocida"
            return "La proteina no se pliega a esta temperatura"
    
    def longitud(self):
        aminoacidos = len(self.aminoacido)
        return f'La proteina tiene {aminoacidos} aminoacidos'
    
    def longitud_secuencia(self):
        total_nucleotidos = len(self.secuencia)
        return f'La proteina tiene {total_nucleotidos} nucleotidos'
    

#probemos 

# gen 
gen1 = Gen("Gen_p68", 100, 250, "ATGGCTAGCTAGCTAGCGCGCGCATATATGCGC")
print(f"Nombre del gen: {gen1.nombre}")
print(f"Secuencia: {gen1.secuencia}")
print(f"Longitud del gen: {gen1.longitud()} bases")
print(f"Contenido GC: {gen1.gc_content():.2%}")

# Secuencia típica de tRNA
secuencia_tRNA = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
trna1 = tRNA("tRNA-Met", 200, 276, secuencia_tRNA, "", "CAU")
print(f"Nombre: {trna1.nombre}")
print(f"Secuencia: {trna1.secuencia}") 
print(f"Anticodón: {trna1.anticodon}")
print(f"Longitud específica de tRNA: {trna1.longitud()}")
print(f"Contenido GC: {trna1.gc_content():.2%}")

# Probar cargar aminoácido
print(trna1.cargar_aminoacido("Metionina"))


# Probar verificación de anticodón
print(f"Anticodón del tRNA: {trna1.anticodon}")

# RNA no codificante
rna_nc = RNA_no_codificante("miRNA-21", 300, 322, "UAGCUUAUCAGACUGAUGUUGA", "silenciamiento génico")
print(f"Nombre: {rna_nc.nombre}")
print(f"Secuencia: {rna_nc.secuencia}")
print(f"Función reguladora: {rna_nc.funcion_reguladora}")
print(f"Longitud: {rna_nc.longitud()} bases")
print(f"Contenido GC: {rna_nc.gc_content():.2%}")
print(f"Estructura: {rna_nc.formar_estructura_secundaria()}")

# Secuencia que codifica para varios aminoácidos
secuencia_proteina = "AUGUUUUCUUAUUGGUGGUGAUAACGCUGCUGCUGAAGCUAAGUAA"
proteina1 = Proteina("Proteína JzWl", 400, 447, secuencia_proteina, "", "UAC", "desconocida")
print(f"Nombre: {proteina1.nombre}")
print(f"Secuencia: {proteina1.secuencia}")
print(f"Anticodón: {proteina1.anticodon}")







