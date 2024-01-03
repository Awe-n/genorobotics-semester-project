from Bio.Blast import NCBIXML
import os
import logging


def parse_blastn_xml(input_xml: str) -> list[tuple[str, int, int, int]] :
    """
    Parses a blastn xml file and returns a list of tuples containing the species name, the coverage percentage and the e-value of the alignment
    """
    with open(input_xml) as file:
        blast_records = NCBIXML.parse(file)
        results = []
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    results.append((alignment.hit_def, hsp.identities/hsp.align_length*100, hsp.expect, hsp.align_length))
        return results
    
def separate_by_species(results: list[tuple[str, int, int, int]]) -> dict[str, list[tuple[int, int, int]]] :
    """
    Separates the results of the blastn into a dictionary with the species name as key and the list of tuples as value
    """
    species = {}
    for result in results:
        species_name = result[0].split()[0] + " " + result[0].split()[1]
        if species_name not in species:
            species[species_name] = []
        species[species_name].append((result[1], result[2], result[3]))
    return species

def reduce_species_results(results: dict[str, list[tuple[int, int, int]]]) -> dict[str, tuple[int, int, int]] :
    """
    Reduces the results of the blastn to the best result for each species
    """
    reduced_results = {}
    for species in results:
        avg_alignment = 0
        avg_evalue = 0
        for result in results[species]:
            avg_alignment += result[0]
            avg_evalue += result[1]
        avg_alignment /= len(results[species])
        avg_evalue /= len(results[species])

        reduced_results[species] = (avg_alignment, avg_evalue)
    return reduced_results


def select_best_species(results: dict[str, tuple[int, int]]) -> tuple[str, tuple[int, int]] :
    """
    Selects the best species from the results of the blastn
    """
    best_species = None
    best_alignment = 0
    best_evalue = float("inf")
    for species in results:
        if results[species][1] < best_evalue:
            best_species = species
            best_alignment = results[species][0]
            best_evalue = results[species][1]
        elif results[species][1] == best_evalue:
            if results[species][0] > best_alignment:
                best_species = species
                best_alignment = results[species][0]
                best_evalue = results[species][1]
    return (best_species, (best_alignment, best_evalue))


def get_best_species_from_xml(input_xml: str) -> tuple[str, tuple[int, int]] :
    """
    Returns the best species from the results of the blastn
    """
    results = parse_blastn_xml(input_xml)
    species = separate_by_species(results)
    reduced_results = reduce_species_results(species)
    return select_best_species(reduced_results)