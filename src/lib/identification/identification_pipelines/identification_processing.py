from Bio.Blast import NCBIXML

def parse_blastn_xml(input_xml: str) -> list[tuple[str, int, int, int]]:
    """
    Parses a blastn XML file and extracts relevant information from it.

    Args:
        input_xml (str): The path to the blastn XML file.

    Returns:
        list[tuple[str, int, int, int]]: A list of tuples containing the hit definition,
        percentage identity, expect value, and alignment length for each alignment in the XML file.
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
    Separates the results by species.

    Args:
        results (list[tuple[str, int, int, int]]): A list of tuples representing the results of a blastn query.
            Each tuple contains the hit definition and three integers representing some data.

    Returns:
        dict[str, list[tuple[int, int, int]]]: A dictionary where the keys are species names
            and the values are lists of tuples representing the data for each species.
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
    Reduce the results of species identification by calculating the average alignment and average e-value for each species.

    Args:
        results (dict[str, list[tuple[int, int, int]]]): A dictionary containing the results of species identification. The keys are species names, and the values are lists of tuples representing alignment, e-value, and alignment length.

    Returns:
        dict[str, tuple[int, int, int]]: A dictionary containing the reduced results of species identification. The keys are species names, and the values are tuples representing the average alignment and average e-value.
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


def select_best_species(results: dict[str, tuple[int, int]]) -> tuple[str, tuple[int, int]]:
    """
    Selects the best species based on the given results (highest alignment, lowest e-value)

    Args:
        results (dict[str, tuple[int, int]]): A dictionary containing species as keys and tuples of alignment and e-value as values.

    Returns:
        tuple[str, tuple[int, int]]: A tuple containing the best species and its corresponding alignment and e-value.
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
    Get the best species from a blastn XML file, by following these steps:
    1. Parse the XML file.
    2. Separate the results by species.
    3. Reduce the results by calculating the average alignment and average e-value for each species.
    4. Select the best species based on the given results (highest alignment, lowest e-value)
    
    Parameters:
    input_xml (str): The path to the input XML file containing blastn results.

    Returns:
    tuple[str, tuple[int, int]]: A tuple containing the best species name and a tuple of two integers representing the start and end positions of the best species in the input XML file.
    """
    results = parse_blastn_xml(input_xml)
    species = separate_by_species(results)
    reduced_results = reduce_species_results(species)
    return select_best_species(reduced_results)