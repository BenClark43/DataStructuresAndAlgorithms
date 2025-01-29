from Bio.Data.IUPACData import protein_letters
from Bio.Align import substitution_matrices
from Bio import SeqIO
import itertools
import ahocorasick
blosum62 = substitution_matrices.load('BLOSUM62')

class Protein:
    def __init__(self, id: str, name: str, record):
        self.id = id
        self.name = name
        self.record = record

class Seed:
    def __init__(self, protein: Protein, kmer: str, loc: int):
        self.protein = protein
        self.kmer = kmer
        self.loc = loc

class Kmer:
    def __init__(self, string, start, end, query, score: int):
        self.string = string
        self.start = start
        self.end = end
        self.query = query
        self.score = score

class Result:
    def __init__(self, protein: Protein, query_start: int, query_end: int, query_text: str, match_start: int, match_end: int, match_text: str, score: int):
        self.protein = protein
        self.query_start = query_start
        self.query_end = query_end
        self.query_text = query_text
        self.match_start = match_start
        self.match_end = match_end
        self.match_text = match_text
        self.score = score

def format_database(database) -> list:
    return [Protein(record.id, record.description, record.seq) for record in database]

def generate_neighbourhood(query: str, k: int, threshold: int) -> dict:
    neighbourhood = {}

    for i in range(len(query) - k + 1):
        kmer = query[i:i+k]
        kmer_score = sum(blosum62[char][char] for char in kmer)
            
        if kmer_score >= threshold: 
            neighbourhood[kmer] = Kmer(kmer, i, i + k - 1, query, kmer_score)

            for mutation in itertools.product(protein_letters, repeat=k):
                mutated_kmer = "".join(mutation)
                mutation_score = sum( blosum62[kmer[idx]] [mutation[idx] ] for idx in range(k))

                if mutation_score >= threshold:
                    if mutated_kmer not in neighbourhood or neighbourhood[mutated_kmer].score < mutation_score:
                        neighbourhood[mutated_kmer] = Kmer(mutated_kmer, i, i+k-1, query, mutation_score)

    return neighbourhood

def seed_search(neighbourhood: dict, database) -> dict:
    seeds = {entry.id: [] for entry in database}

    if neighbourhood:

        automaton = ahocorasick.Automaton()
        for kmer in neighbourhood.keys():
            automaton.add_word(kmer, kmer)
        automaton.make_automaton()

        for entry in database:
            matches = {}
            for end_index, pattern in automaton.iter(str(entry.record)):
                if pattern not in matches:
                    matches[pattern] = []
                matches[pattern].append(end_index)

            seeds[entry.id] = [
                Seed(entry, neighbourhood[pattern], location)
                for pattern, locations in matches.items()
                for location in locations
            ]

    return {entry_id: seed_list for entry_id, seed_list in seeds.items() if seed_list}

def extend_seed(seed: Seed, minimum_score: int, query: str, text: str) -> Result:
    score = seed.kmer.score
    kmer_start, kmer_end = seed.kmer.start, seed.kmer.end
    match_start, match_end = seed.loc - len(seed.kmer.string)+1, seed.loc
    query_len, text_len = len(query), len(text)

    can_extend_start = lambda: kmer_start > 0 and match_start > 0
    can_extend_end = lambda: kmer_end < query_len - 1 and match_end < text_len - 1

    while (score >= minimum_score and (can_extend_start() or can_extend_end())):
        if can_extend_start():
            kmer_start -= 1
            match_start -= 1
            score += blosum62[query[kmer_start]][text[match_start]]
        if can_extend_end():
            kmer_end += 1
            match_end += 1
            score += blosum62[query[kmer_end]][text[match_end]]

    return Result(
        protein = seed.protein.name,
        query_start = kmer_start,
        query_end = kmer_end,
        query_text = query[kmer_start:kmer_end+1],
        match_start = match_start,
        match_end = match_end,
        match_text = text[match_start:match_end+1],
        score = score
    )

def identify_alignments(seed_dict: dict, minimum_score: int, threshold: int) -> dict:
    extended_seeds = {}
    high_scoring_alignments = []
    
    for seed_list in seed_dict.values():
        extended_seeds[seed_list[0].protein.id] = []
        hsas = []
 
        for seed in seed_list:
            alignment = extend_seed(seed, minimum_score, seed.kmer.query, seed.protein.record)
            
            if alignment.score >= threshold:
                if all(len(set(range(alignment.match_start, alignment.match_end)).intersection(range(value.match_start, value.match_end))) == 0 for value in hsas):
                    hsas.append(alignment)  
        high_scoring_alignments.extend(hsas)

    return high_scoring_alignments

"""
    database: a fasta file containing the proteins
    query: the target string to search for (must be protein characters)
    kmer: the size of the kmers in the neighbourhood
    neighbourhood: the threshold that additional kmers need to pass to be added into the neighbourhood (should be set relative to kmer)
    extension: the amount the score of a match can fall before it stops being extended further
    minimum_alignment_score : The minimum score a match needs to pass in order to be returned
"""
def blast(database: str, query: str, kmer: int, neighbourhood: int, extension: int, minimum_alignment_score: int):
    try:
        database_iterator = SeqIO.parse(database, "fasta")
        proteins =  [Protein(record.id, record.description, record.seq) for record in database_iterator]

        kmers = generate_neighbourhood(query, kmer, neighbourhood)
        seeds = seed_search(kmers, proteins)
        high_scoring_alignments = identify_alignments(seeds, neighbourhood - extension, minimum_alignment_score)
        high_scoring_alignments.sort(key=lambda alignment: alignment.score, reverse=True)

        print(" ---- High Scoring Alignments ---- ")
        for alignment in high_scoring_alignments:
            print("{}".format(alignment.protein))
            print(" * Score: {}".format(alignment.score))
            print(" * Record sequence: {} [{}, {}]".format(alignment.match_text, alignment.match_start, alignment.match_end))
            print(" * Query  sequence: {} [{}, {}]".format(alignment.query_text, alignment.query_start, alignment.query_end))   


    except FileNotFoundError:
        print("The file does not exist")

def example():
    blast('proteins.fasta', 'SLKGGKAHWYVPDANKLRSA', 3, 14, 5, 36)

"""
    ---- High Scoring Alignments ----
    Protein_Sequence_5 | hypothetical protein B
    * Score: 109.0
    * Record sequence: SLKGGKAHWYVPDANKLRSA [324, 343]
    * Query  sequence: SLKGGKAHWYVPDANKLRSA [0, 19]
    Protein_Sequence_2 | hypothetical protein Y
    * Score: 98.0
    * Record sequence: SLKTGKAHWYVPDANKLRSS [324, 343]
    * Query  sequence: SLKGGKAHWYVPDANKLRSA [0, 19]
    Protein_Sequence_3 | hypothetical protein Z
    * Score: 94.0
    * Record sequence: SLKTGRAHWYVPDTNKLRSA [324, 343]
    * Query  sequence: SLKGGKAHWYVPDANKLRSA [0, 19]
    Protein_Sequence_4 | hypothetical protein A
    * Score: 93.0
    * Record sequence: SLKTGKAHWYVPDTNKLRST [324, 343]
    * Query  sequence: SLKGGKAHWYVPDANKLRSA [0, 19]
    Protein_Sequence_1 | hypothetical protein X
    * Score: 92.0
    * Record sequence: SLKFGKTHWYVPDANKLRST [324, 343]
    * Query  sequence: SLKGGKAHWYVPDANKLRSA [0, 19]
"""

if __name__ == '__main__':
    example()