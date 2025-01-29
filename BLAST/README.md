# Basic Local Alignment Search Tool (BLAST),

##### Set up
Requiremented libraries can be installed using `pip`
```
pip install -r requirements.txt
```

#### To run program
```
python ./blast.py
```
This will run the example query with the default settings on the example datasource

## Algorithm
##### 1. Blosum Score
BLAST works with patterns that do not exactly match the text. These mismatches arise from mutations in proteins. Since some mutations are more common than others, the algorithm uses a BLOSUM score to assign a weight to each mutation based on its likelihood. These scores are then used to determine which missmatches are better than others.

##### 2. Neighbourhood Generation
In order to make the query string more mangageable, it is devided into smaller **k-mers** (substrings of length K).
For example, with `K=5`:
`ABCDEFGH` -> `ABCDE`, `BCDEF`, `CDEFG`, `DEFGH`

In practise, when `K` is large we want to allow for small mismatches to generate more seeds.
This is achieved by creating a **neighborhood** of similar k-mers. These are derived by mutating the input k-mers, prioritizing mutations with high BLOSUM scores.

##### 3. Seed Search
After the neighbourhood of kmers can been constructed, the entire text is searched. This step is done using an **exact pattern matching** algorithm (Aho Corasick). 
The resulting matches are called seeds

##### 4. Seed-Extension and High-Scoring Alignments
Each seed keeps track of its current match score (determined by blosum values) which will initial be quite high. The seeds are then extended in each direction until the the match score drops below a threshold.

After all seeds have been extended, large seeds with high match scores are determined to be high values alignments, and are returned to the user.

<br>

## Aho Corasick (Exact Pattern Matching)
#### Algorithm 
Aho Corasick allows for multiple patterns to be search for simulationaly in a single pass through the text. 

This is accomplished by building a Trie. Each node contains a single letter from an input. Then, as symbols are read from the text, a pointer moves through the trie to determine if a match has been achieved. 

The Trie contains additional pointers between the nodes reffered to as arcs which support multiple matches occuring from a single character read.
<br>
<img src="https://upload.wikimedia.org/wikipedia/commons/9/90/A_diagram_of_the_Aho-Corasick_string_search_algorithm.svg" alt="Aho-Corasick Diagram" width="340">

