#Genere la clase gen con al menos 3 atributos y 3 métodos (contando el contructor)

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

#Provemos
DNA_sequence = Gen("luxi", 10, 18, "ATTTGCGA")
print(f'La longitud total considerando su posicion genomica es de: {DNA_sequence.longitud()}')
print(f'El porcentaje de GC de la secuencia es del: {DNA_sequence.gc_content()*100}%')