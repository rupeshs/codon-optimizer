import csv
import python_codon_tables as pct
from dnachisel import DnaOptimizationProblem, EnforceTranslation, EnforceGCContent, CodonOptimize
from prettytable import PrettyTable
import numpy

# Setting seed, DNA Chisel Reproducibility
numpy.random.seed(10110)

class VaccineCodonOptimiser():
    def __init__(self):
        self.__virus_codons = []
        self.__vaccine_codons = []
        self.__vaccine_codons_gen = []
        self.__codon_match_ratio = 0
        self.__neucleotide_match_ratio = 0
        self.__amino_acids = {}
        self.__load_amnioacid_mappings()

    def load_codons(self, side_by_side_file):
        with open(side_by_side_file) as csvfile:
            data = csv.DictReader(csvfile)
            for row in data:
                self.__virus_codons.append(row['codonOrig'])
                self.__vaccine_codons.append(row['codonVaccine'])

    def __load_amnioacid_mappings(self):
        with open("codon-table-grouped.csv") as csvfile:
            data = csv.DictReader(csvfile)
            for row in data:
                self.__amino_acids[row['codon']] = row["aminoacid"]

    def get_vaccine_codons(self):
        return self.__vaccine_codons_gen

    def get_virus_codons(self):
        return self.__virus_codons

    def get_codon_match_ratio(self):
        match_count = 0
        match_count = sum(x == y for x, y in zip(
            self.__vaccine_codons_gen, self.__vaccine_codons))
        self.__codon_match_ratio = (
            match_count / len(self.__vaccine_codons)) * 100
        return round(self.__codon_match_ratio, 2)

    def get_nucleotide_match_ratio(self):
        match = 0
        total_len = 0
        for x, y in zip(self.__vaccine_codons, self.__vaccine_codons_gen):
            for na, nb in zip(x, y):
                if na == nb:
                    match += 1
            total_len += 3
        ratio = (match / total_len) * 100.0
        return round(ratio, 2)

    def get_gc_ratio(self):
        #GC-content percentage calculated as per https://en.wikipedia.org/wiki/GC-content
        g_count = 0
        c_count = 0
        total_nt = 0
        for codon in self.__vaccine_codons_gen:
            total_nt += 1
            for nt in codon:
                if nt == "G":
                    g_count += 1
                elif nt == "C":
                    c_count += 1
        gc_ratio = (g_count + c_count) / (total_nt * 3) * 100
        return round(gc_ratio, 2)

    def optimize_frequent(self, codon_table):
        table = pct.get_codons_table(codon_table)
        self.__vaccine_codons_gen.clear()
        for codon in self.__virus_codons:
            amino_acid = self.__amino_acids[codon]
            if amino_acid == "s":
                amino_acid = "*"
            sorted_amino = dict(
                sorted(table[amino_acid].items(), key=lambda item: item[1], reverse=True))
            new_codon = list(sorted_amino.keys())[0]
            self.__vaccine_codons_gen.append(new_codon)
        return

    def get_strand(self, codons):
        codon_strand = ""
        for codon in codons:
            codon_strand = codon_strand + codon
        return codon_strand

    def optimize(self, codon_table):
        self.optimize_frequent(codon_table)
        # return
        opt_codons = self.__vaccine_codons_gen.copy()
        self.__vaccine_codons_gen.clear()
        vac_strand = self.get_strand(opt_codons)
        #vir_strand = self.get_strand(self.__virus_codons)
        codon_table = pct.get_codons_table(codon_table)
        problem = DnaOptimizationProblem(
            sequence=vac_strand,
            constraints=[
                EnforceTranslation(genetic_table='Standard',
                                   start_codon='ATG'),
                EnforceGCContent(mini=0.54, maxi=0.9, window=120)
            ],
            objectives=[
                CodonOptimize(method="use_best_codon",
                              codon_usage_table=codon_table)
            ]
        )
        problem.resolve_constraints()
        problem.optimize()
        self.__vaccine_codons_gen = []
        count = 1
        vcodon = ""
        for x in problem.sequence:
            if count % 3 == 0:
                vcodon += x
                self.__vaccine_codons_gen.append(vcodon)
                vcodon = ""
            else:
                vcodon += x
            count += 1
        return
    
if __name__ == "__main__":
    
    vacc_optimizer = VaccineCodonOptimiser()
    vacc_optimizer.load_codons("side-by-side.csv")
    
    ptbl = PrettyTable()
    ptbl.field_names = ["Species", "Codon Match %",
                        "Nucleotide Match %", "GC ratio %"]
    species = ["h_sapiens_9606", "m_musculus_10090"]
    for speci in species:
        vacc_optimizer.optimize(speci)
        ptbl.add_row([speci,
                      vacc_optimizer.get_codon_match_ratio(),
                      vacc_optimizer.get_nucleotide_match_ratio(),
                      vacc_optimizer.get_gc_ratio()])
    print(ptbl)
    #print(vacc_optimizer.get_strand(vacc_optimizer.get_vaccine_codons()))
