from neo4j import GraphDatabase
from collections import deque
from rdkit import Chem


URI = "bolt://localhost:7687"
USER = "neo4j"
PASSWORD = "chem_react"

# Retrosynthesis via BFS with batched database queries
# ------------------------------

def expand_precursors_batch(target_smiles_list, driver):
    """
    Finds all reactions that produce any of the target molecules in the list.
    Returns a dictionary mapping each target to its list of precursor reactions.
    """
    query = """
    MATCH (precursor:Molecule)-[:REACTANT_IN]->(r:Reaction)-[:PRODUCES]->(m:Molecule)
    WHERE m.smiles IN $targets
    RETURN m.smiles AS target, 
           r.reaction_id AS rid,
           collect(DISTINCT precursor.smiles) AS precursors
    """
    
    precursor_map = {target: [] for target in target_smiles_list}
    
    with driver.session() as session:
        results = session.run(query, targets=target_smiles_list)
        for record in results:
            target = record["target"]
            precursor_map[target].append((record["rid"], record["precursors"]))
            
    return precursor_map



def retrosynthesis_bfs(target, building_blocks, driver, max_depth=10):
    """
    Performs a Breadth-First Search for a retrosynthetic path using
    level-wise expansion and batched database queries.
    Returns a list of reaction IDs representing a valid path.
    """
    # Get the element sets for building blocks and target
    building_blocks_elems = set()
    for bb in building_blocks:
        building_blocks_elems.update(get_elems(bb))
    target_elems = get_elems(target)

    # First feasibility check based on elemental composition
    if not contains_elements(building_blocks_elems, target_elems):
        print("No reaction path possible, because some elements in the target are not present in the building blocks.")
        return None
    
    # Queue stores tuples of (set_of_molecules_to_synthesize, path_of_reaction_ids)
    queue = deque([ (frozenset([target]), []) ])
    # Visited set stores frozensets of molecules to avoid redundant work and cycles
    visited = {frozenset([target])}

    for depth in range(max_depth):
        # First, expand to precursors for all molecules at this level in a batch
        mols_to_expand_this_level = set()
        level_size = len(queue)
        
        if level_size == 0:
            break # No more paths to explore

        # Use a temporary list to hold states for the current level processing batch
        current_level_states = [queue.popleft() for _ in range(level_size)]

        for current_molecules, _ in current_level_states:
            # Check if this entire state is viable
            # If any molecule has "impossible" elements, this branch is dead.
            is_state_viable = True
            for m in current_molecules:
                if m not in building_blocks and not get_elems(m).issubset(building_blocks_elems):
                    is_state_viable = False
                    break

            if not is_state_viable:
                continue # Skip this state entirely, don't expand any molecules in it

            # If viable, proceed to collect molecules for batching...
            mols_to_expand_this_level.update(
                m for m in current_molecules 
                if m not in building_blocks) # We already checked subset condition above

        # --- Batch Query ---
        # If there are no molecules to expand on this level, we might have found our solution
        if not mols_to_expand_this_level:
            # Check if any of the current states are solved
            for current_molecules, path in current_level_states:
                if current_molecules.issubset(building_blocks):
                    return path # Should have been caught earlier, but as a safeguard
            continue
            
        precursor_map = expand_precursors_batch(list(mols_to_expand_this_level), driver)

        # --- Level Expansion ---
        # Now process each state from the current level
        for current_molecules, path in current_level_states:

            # Try to expand each molecule in the current set
            for mol_to_expand in current_molecules:
                if mol_to_expand in building_blocks:
                    continue

                # Use the pre-fetched results from our batch query
                expansions = precursor_map.get(mol_to_expand, [])

                for reaction_id, precursors in expansions:
                    # Create the new set of molecules for the next state
                    next_molecules = (current_molecules - {mol_to_expand}).union(precursors)
                    
                    # Check if all molecules are building blocks (goal condition)
                    if next_molecules.issubset(building_blocks):
                        # A complete path is found
                        return path + [reaction_id]

                    if next_molecules not in visited:
                        visited.add(next_molecules)
                        new_path = path + [reaction_id]
                        queue.append((next_molecules, new_path))
    return None

# Hard constraints for feasibility checkng to prune search space
# ------------------------------
# Simple chemical constraint: elements cannot be created out of nothing
# -> Elements(target) ⊆ Elements(reactants)

def get_elems(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        return set() # Or raise error / handle gracefully
    return {a.GetSymbol() for a in mol.GetAtoms()}

def contains_elements(reactant_elems, product_elems):
    return product_elems.issubset(reactant_elems)

# Getting the full synthesis subgraph
# ------------------------------
def get_synthesis_subgraph(reaction_ids, driver):
    """
    Given a list of reaction IDs, returns the full subgraph containing
    all molecules and reactions in the synthesis pathway.
    """
    # This query uses pattern comprehensions, which is a more modern and
    # efficient syntax that avoids the deprecated CALL subquery.
    query = """
    MATCH (r:Reaction)
    WHERE r.reaction_id IN $rids
    RETURN r.reaction_id AS reaction,
           [(reactant:Molecule)-[:REACTANT_IN]->(r) | reactant.smiles] AS reactants,
           [(r)-[:PRODUCES]->(product:Molecule) | product.smiles] AS products
    """
    with driver.session() as session:
        results = session.run(query, rids=reaction_ids)
        return [dict(record) for record in results]


def get_reaction_procedures(reaction_ids, driver):
    """
    Given a list of reaction IDs (ordered), returns a list of dictionaries
    containing the reaction_id and its procedure text.
    The output list preserves the order of the input reaction_ids.
    """
    query = """
    MATCH (r:Reaction)
    WHERE r.reaction_id = $rid
    RETURN r.reaction_id AS reaction_id, r.procedure AS procedure
    """
    
    procedures = []
    with driver.session() as session:
        for rid in reaction_ids:
            result = session.run(query, rid=rid).single()
            if result:
                procedures.append(dict(result))
            else:
                procedures.append({"reaction_id": rid, "procedure": "No procedure details found."})
                
    return procedures

# ------------------------------

def main(target_molecule, starting_molecules):

    driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))

    print(f"Searching for a synthesis path for: {target_molecule}")
    
    # Find a path (list of reaction_ids)
    path_rids = retrosynthesis_bfs(target_molecule, starting_molecules, driver)

    if path_rids:
        print(f"Found a path with {len(path_rids)} reaction(s):")
        print(path_rids)
        
        # Retrieve the full subgraph for the path
        subgraph = get_synthesis_subgraph(path_rids, driver)
        
        print("\nSynthesis Subgraph:")
        for reaction_step in subgraph:
            print(f"  Reaction: {reaction_step['reaction']}")
            print(f"    Reactants: {reaction_step['reactants']}")
            print(f"    Products: {reaction_step['products']}")
    else:
        print("No synthesis path found within the maximum depth.")

    driver.close()


if __name__ == "__main__":
    target_molecule = "COc1cc2c(cc1Nc1ncc(C)c(-c3cnn4ccccc34)n1)N(C(C)=O)CC2"
    starting_molecules = set([
        "Cc1cnc(Cl)nc1-c1cnn2ccccc12",
        "CCO",
        "COc1cc2c(cc1[N+](=O)[O-])N(C(C)=O)CC2"
    ])

    main(target_molecule, starting_molecules)
